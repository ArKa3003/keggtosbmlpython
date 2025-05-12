#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import requests
import argparse
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set, Union
import xml.etree.ElementTree as ET
from xml.dom import minidom
import time
import logging

# --- Basic Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class KEGGtoSBMLConverter:
    """
    Convert KEGG pathway data to SBML format.
    Fetches pathway, compound, and reaction data from the KEGG REST API.
    """

    def __init__(self, pathway_id: str = "map00790", delay: float = 0.1):
        """
        Initialize the converter with a KEGG pathway ID.

        Args:
            pathway_id: KEGG pathway ID (e.g., map00790 - Folate Biosynthesis)
            delay: Delay in seconds between KEGG API calls to avoid overloading the server.
        """
        if not pathway_id.startswith("map") and not pathway_id.startswith("ko"):
             # Assume it might be a search term or specific pathway code
             logging.warning(f"Pathway ID '{pathway_id}' might not be a standard map ID. Proceeding anyway.")
        # Basic validation - ensure it's not empty
        if not pathway_id:
            raise ValueError("Pathway ID cannot be empty.")

        self.pathway_id = pathway_id
        self.base_url = "http://rest.kegg.jp"
        self.api_delay = delay # seconds between API calls
        self.species = {}  # Dictionary to store species/compounds: {kegg_id: data_dict}
        self.reactions = {}  # Dictionary to store reactions: {kegg_id: data_dict}
        self.pathway_name = ""
        self.compartments = {"default": {"name": "cytosol", "size": 1.0}} # Default compartment

        # For managing unique SBML IDs
        self.species_sbml_ids = set()
        self.reaction_sbml_ids = set()
        self.parameter_sbml_ids = set()
        self.compartment_sbml_ids = set()

    def _make_api_request(self, url: str) -> Optional[str]:
        """
        Makes a request to the KEGG API with error handling and delay.

        Args:
            url: The URL to fetch.

        Returns:
            The response text if successful, None otherwise.
        """
        try:
            time.sleep(self.api_delay) # Respect KEGG API usage policies
            response = requests.get(url, timeout=30) # Added timeout
            response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
            return response.text
        except requests.exceptions.RequestException as e:
            logging.error(f"API request failed for {url}: {e}")
            return None
        except Exception as e:
            logging.error(f"An unexpected error occurred during API request for {url}: {e}")
            return None

    def _generate_unique_sbml_id(self, base_id: str, existing_ids: Set[str]) -> str:
        """
        Generates a unique, SBML-compliant ID.

        Args:
            base_id: The desired base for the ID (will be cleaned).
            existing_ids: A set of already used IDs in the current scope.

        Returns:
            A unique SBML ID.
        """
        # SBML IDs must start with a letter or underscore, and contain only letters, numbers, underscores
        # Replace invalid characters with underscore
        sbml_base = re.sub(r'^[^a-zA-Z_]', '_', base_id)
        sbml_base = re.sub(r'[^a-zA-Z0-9_]', '_', sbml_base)

        if not sbml_base: # Handle cases where the base_id was only invalid chars
            sbml_base = "_id"

        unique_id = sbml_base
        counter = 1
        while unique_id in existing_ids:
            unique_id = f"{sbml_base}_{counter}"
            counter += 1
        existing_ids.add(unique_id)
        return unique_id

    def fetch_pathway_data(self) -> Optional[str]:
        """
        Fetch main pathway data (structure, components) from KEGG API.

        Returns:
            Raw KEGG pathway data as string, or None on failure.
        """
        logging.info(f"Fetching pathway data for {self.pathway_id}...")
        url = f"{self.base_url}/get/{self.pathway_id}"
        data = self._make_api_request(url)
        if data:
            logging.info(f"Successfully fetched pathway data for {self.pathway_id}.")
        else:
            logging.error(f"Failed to fetch pathway data for {self.pathway_id}.")
        return data

    def fetch_pathway_name(self, pathway_data: Optional[str]) -> str:
        """
        Extract the name of the pathway from the raw pathway data or fetch it.

        Args:
            pathway_data: Raw pathway data string (optional).

        Returns:
            Pathway name, or a default name if not found.
        """
        if pathway_data:
            # Try extracting from the NAME field in the main pathway data first
            name_match = re.search(r"^NAME\s+(.*)", pathway_data, re.MULTILINE)
            if name_match:
                self.pathway_name = name_match.group(1).strip()
                logging.info(f"Found pathway name from data: {self.pathway_name}")
                return self.pathway_name

        # Fallback to separate API call if not found in main data or data is missing
        logging.info(f"Fetching pathway name separately for {self.pathway_id}...")
        url = f"{self.base_url}/list/pathway" # More robust way to get all names
        data = self._make_api_request(url)
        if data:
            for line in data.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2 and parts[0] == f"path:{self.pathway_id}":
                    # Example: path:map00790	Folate biosynthesis - Homo sapiens (human)
                    # Take the part before species information if present
                    name = parts[1].split(' - ')[0].strip()
                    self.pathway_name = name
                    logging.info(f"Found pathway name via API list: {self.pathway_name}")
                    return self.pathway_name

        logging.warning(f"Could not determine pathway name for {self.pathway_id}. Using default.")
        self.pathway_name = f"KEGG_{self.pathway_id}"
        return self.pathway_name

    def parse_compound_data(self, compound_id: str) -> Optional[Dict]:
        """
        Fetch and parse compound data. Caches results in self.species.

        Args:
            compound_id: KEGG compound ID (e.g., C00001).

        Returns:
            Dictionary with compound information, or None on failure.
        """
        if not compound_id or not compound_id.startswith('C'):
            logging.warning(f"Invalid compound ID format: {compound_id}. Skipping.")
            return None

        # Return cached data if available
        if compound_id in self.species:
            return self.species[compound_id]

        logging.info(f"Fetching data for compound {compound_id}...")
        url = f"{self.base_url}/get/{compound_id}"
        data = self._make_api_request(url)
        if not data:
            logging.error(f"Failed to fetch data for compound {compound_id}.")
            # Store minimal info to avoid re-fetching repeatedly
            self.species[compound_id] = {"id": compound_id, "name": compound_id, "kegg_id": compound_id}
            return self.species[compound_id]

        compound_info = {"id": compound_id, "kegg_id": compound_id, "name": compound_id, "formula": None, "charge": None, "names": []}

        # Use more robust regex for potentially multi-line fields ending before the next keyword
        name_match = re.search(r"^NAME\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if name_match:
            # Get all names, primary is usually the first line
            raw_names = name_match.group(1).strip()
            names_list = [n.strip().rstrip(';') for n in raw_names.split('\n')]
            compound_info["name"] = names_list[0] if names_list else compound_id
            compound_info["names"] = names_list

        formula_match = re.search(r"^FORMULA\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if formula_match:
            compound_info["formula"] = formula_match.group(1).strip()

        charge_match = re.search(r"^CHARGE\s+([-+]?\d+)", data, re.MULTILINE)
        if charge_match:
             try:
                compound_info["charge"] = int(charge_match.group(1))
             except ValueError:
                logging.warning(f"Could not parse charge for {compound_id} from: {charge_match.group(1)}")


        logging.debug(f"Parsed compound {compound_id}: {compound_info['name']}, Formula: {compound_info['formula']}")
        self.species[compound_id] = compound_info
        return compound_info


    def parse_reaction_data(self, reaction_id: str) -> Optional[Dict]:
        """
        Fetch and parse reaction data. Caches results in self.reactions.

        Args:
            reaction_id: KEGG reaction ID (e.g., R00001).

        Returns:
            Dictionary with reaction information, or None on failure.
        """
        if not reaction_id or not reaction_id.startswith('R'):
            logging.warning(f"Invalid reaction ID format: {reaction_id}. Skipping.")
            return None

        # Return cached data if available
        if reaction_id in self.reactions:
             return self.reactions[reaction_id]

        logging.info(f"Fetching data for reaction {reaction_id}...")
        url = f"{self.base_url}/get/{reaction_id}"
        data = self._make_api_request(url)
        if not data:
            logging.error(f"Failed to fetch data for reaction {reaction_id}.")
             # Store minimal info
            self.reactions[reaction_id] = {"id": reaction_id, "name": reaction_id, "kegg_id": reaction_id}
            return self.reactions[reaction_id]


        reaction_info = {
            "id": reaction_id,
            "kegg_id": reaction_id,
            "name": reaction_id,
            "equation_str": None,
            "definition_str": None, # Raw DEFINITION field might be more consistent
            "substrates": [], # List of tuples (compound_kegg_id, stoichiometry)
            "products": [],   # List of tuples (compound_kegg_id, stoichiometry)
            "enzymes": [],    # List of KEGG enzyme IDs (e.g., EC numbers)
            "reversible": False
        }

        # Extract name
        name_match = re.search(r"^NAME\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if name_match:
            reaction_info["name"] = name_match.group(1).strip().split('\n')[0].strip('; ') # Often has trailing ';'

        # Extract definition (often cleaner than EQUATION)
        definition_match = re.search(r"^DEFINITION\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if definition_match:
             reaction_info["definition_str"] = definition_match.group(1).strip()

        # Extract equation (as fallback or for display)
        equation_match = re.search(r"^EQUATION\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if equation_match:
             reaction_info["equation_str"] = equation_match.group(1).strip()

        # Determine reversibility (from DEFINITION or EQUATION)
        eq_str = reaction_info["definition_str"] or reaction_info["equation_str"]
        if eq_str:
            reaction_info["reversible"] = "<=>" in eq_str

        # Extract Enzymes
        enzyme_match = re.search(r"^ENZYME\s+(.*?)(?=\n[A-Z][A-Z]+\s|\Z)", data, re.DOTALL | re.MULTILINE)
        if enzyme_match:
            enzyme_text = enzyme_match.group(1).strip()
            # Match EC numbers (e.g., 1.1.1.1) or KO numbers (K#####) if needed
            reaction_info["enzymes"] = re.findall(r'\d+\.\d+\.\d+\.\d+|\bK\d{5}\b', enzyme_text)


        # --- Parse Substrates and Products from DEFINITION or EQUATION ---
        # Prefer DEFINITION as it's generally cleaner
        target_eq = reaction_info["definition_str"] or reaction_info["equation_str"]

        if target_eq:
            logging.debug(f"Parsing equation for {reaction_id}: {target_eq}")
            try:
                # Split into sides
                if "<=>" in target_eq:
                    parts = target_eq.split("<=>")
                    if len(parts) == 2:
                       substrates_str, products_str = parts
                    else: # Malformed? Log and skip parsing stoichiometry
                        logging.warning(f"Reaction {reaction_id}: Unexpected number of parts with <=>: {target_eq}")
                        substrates_str, products_str = "", ""
                elif "=>" in target_eq:
                    parts = target_eq.split("=>")
                    if len(parts) == 2:
                       substrates_str, products_str = parts
                    else: # Malformed?
                       logging.warning(f"Reaction {reaction_id}: Unexpected number of parts with =>: {target_eq}")
                       substrates_str, products_str = "", ""
                elif "<=" in target_eq: # Handle reverse arrows sometimes seen
                    parts = target_eq.split("<=")
                    if len(parts) == 2:
                       products_str, substrates_str = parts # Swap sides
                       logging.warning(f"Reaction {reaction_id}: Used reverse arrow <=, swapping substrate/product: {target_eq}")
                    else:
                       logging.warning(f"Reaction {reaction_id}: Unexpected number of parts with <=: {target_eq}")
                       substrates_str, products_str = "", ""
                else:
                    # Assume it might be just products if no arrow found (less common)
                    logging.warning(f"Reaction {reaction_id}: No standard arrow found in equation: {target_eq}. Assuming all products.")
                    substrates_str = ""
                    products_str = target_eq

                # Regex to find "([stoichiometry]) [CompoundID]" pairs (stoichiometry optional)
                # Example: "2 C00001 + C00002" -> find (2, C00001), (1, C00002)
                participant_regex = re.compile(r"(?:^|\s*\+\s*)(?:(\d+)\s+)?(C\d{5})\b")

                # Process substrates
                for match in participant_regex.finditer(substrates_str.strip()):
                    stoich_str, comp_id = match.groups()
                    stoich = int(stoich_str) if stoich_str else 1
                    reaction_info["substrates"].append((comp_id, stoich))
                    # Ensure compound data is fetched later if not already present
                    self.parse_compound_data(comp_id)


                # Process products
                for match in participant_regex.finditer(products_str.strip()):
                    stoich_str, comp_id = match.groups()
                    stoich = int(stoich_str) if stoich_str else 1
                    reaction_info["products"].append((comp_id, stoich))
                     # Ensure compound data is fetched later if not already present
                    self.parse_compound_data(comp_id)


                logging.debug(f"Parsed Reaction {reaction_id}: Subs={reaction_info['substrates']}, Prods={reaction_info['products']}, Rev={reaction_info['reversible']}")

            except Exception as e:
                logging.error(f"Error parsing equation for reaction {reaction_id}: {target_eq} - {e}")

        else:
             logging.warning(f"Reaction {reaction_id}: Could not find DEFINITION or EQUATION field.")

        self.reactions[reaction_id] = reaction_info
        return reaction_info

    def parse_pathway(self) -> bool:
        """
        Parse the main KEGG pathway data file to identify components.
        Fetches detailed data for each component.

        Returns:
            True if parsing was successful (found some reactions/species), False otherwise.
        """
        pathway_data = self.fetch_pathway_data()
        if not pathway_data:
            return False # Stop if we can't get the main file

        self.fetch_pathway_name(pathway_data) # Try extracting name from file first

        # Process Compounds listed explicitly in the pathway file
        # These are often just for visual map layout, reactions are the source of truth
        # Using robust section extraction looking for next keyword or end of file
        compound_pattern = r"^COMPOUND\s+(.*?)(?=\n^[A-Z][A-Z]+\s|\Z)"
        compound_match = re.search(compound_pattern, pathway_data, re.DOTALL | re.MULTILINE)
        compounds_in_pathway = set()
        if compound_match:
            compounds_text = compound_match.group(1).strip()
            for line in compounds_text.split('\n'):
                parts = line.strip().split()
                if parts and parts[0].startswith('C'):
                    compound_id = parts[0]
                    compounds_in_pathway.add(compound_id)
                    if compound_id not in self.species: # Avoid redundant API calls if reaction parsing already found it
                         self.parse_compound_data(compound_id)
            logging.info(f"Found {len(compounds_in_pathway)} compounds explicitly listed in pathway map.")
        else:
             logging.warning(f"Could not find COMPOUND section in {self.pathway_id} data.")


        # Process Reactions listed explicitly in the pathway file
        # This section links Reaction IDs to pathway context (like enzymes)
        reaction_pattern = r"^REACTION\s+(.*?)(?=\n^[A-Z][A-Z]+\s|\Z)"
        reaction_match = re.search(reaction_pattern, pathway_data, re.DOTALL | re.MULTILINE)
        reactions_in_pathway = set()
        if reaction_match:
            reactions_text = reaction_match.group(1).strip()
            reaction_ids_from_pathway = set(re.findall(r'\b(R\d{5})\b', reactions_text))

            logging.info(f"Found {len(reaction_ids_from_pathway)} reaction IDs in pathway REACTION section. Fetching details...")
            for reaction_id in reaction_ids_from_pathway:
                reactions_in_pathway.add(reaction_id)
                reaction_data = self.parse_reaction_data(reaction_id)
                # Ensure species from reactions are parsed (parse_reaction_data handles this now)
                # if reaction_data:
                #    for comp_id, _ in reaction_data.get("substrates", []) + reaction_data.get("products", []):
                #        if comp_id not in self.species:
                #            self.parse_compound_data(comp_id)
            logging.info(f"Finished processing details for {len(reactions_in_pathway)} reactions.")

        else:
            logging.warning(f"Could not find REACTION section in {self.pathway_id} data.")

        # Simple check: if we didn't find any species or reactions, something likely went wrong
        if not self.species or not self.reactions:
            logging.error("Parsing finished, but no species or reactions were successfully parsed. Check pathway ID and KEGG format.")
            return False

        logging.info(f"Pathway parsing complete. Found {len(self.species)} unique compounds and {len(self.reactions)} unique reactions.")
        return True


    def create_sbml_model(self) -> ET.Element:
        """
        Create SBML model (as XML ElementTree) from parsed KEGG data.

        Returns:
            XML Element containing the SBML model.
        """
        logging.info("Creating SBML model structure...")
        # Create SBML root element
        sbml_attrib = {
            "xmlns": "http://www.sbml.org/sbml/level3/version2/core", # Use L3V2 commonly supported
            "level": "3",
            "version": "2"
        }
        sbml = ET.Element("sbml", sbml_attrib)

        # Create model element
        model_id = self._generate_unique_sbml_id(f"model_{self.pathway_id}", set()) # Model ID needs uniqueness too potentially
        model = ET.SubElement(sbml, "model", {"id": model_id, "name": self.pathway_name})

        # Add notes about the model creation
        notes_str = f"""
        <notes>
          <html xmlns="http://www.w3.org/1999/xhtml">
            <head>
               <title>Model Notes</title>
            </head>
            <body>
              <p>This model was automatically generated from KEGG pathway <a href="https://www.kegg.jp/dbget-bin/www_bget?{self.pathway_id}">{self.pathway_id}</a> ({self.pathway_name}) using KEGGtoSBMLConverter.</p>
              <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
              <p>KEGG API Base URL: {self.base_url}</p>
              <p><b>Disclaimer:</b> This model represents the stoichiometry and connectivity found in KEGG. Kinetic laws are generic placeholders (mass action) and may need refinement for simulation.</p>
             </body>
          </html>
        </notes>
        """
        # Append notes as parsed XML
        try:
            notes_element = ET.fromstring(notes_str)
            model.append(notes_element)
        except ET.ParseError as e:
            logging.warning(f"Could not parse notes XML: {e}")
            # Add as simple text comment instead
            model.append(ET.Comment(f"Model generated from KEGG {self.pathway_id} ({self.pathway_name})"))

        # --- Create Compartments ---
        listOfCompartments = ET.SubElement(model, "listOfCompartments")
        for comp_kegg_id, comp_data in self.compartments.items():
            comp_sbml_id = self._generate_unique_sbml_id(comp_kegg_id, self.compartment_sbml_ids)
            comp_data['sbml_id'] = comp_sbml_id # Store for reference
            compartment = ET.SubElement(listOfCompartments, "compartment", {
                "id": comp_sbml_id,
                "name": comp_data.get("name", comp_kegg_id),
                "constant": "true", # Assume constant volume
                "spatialDimensions": "3", # Assume 3D
                "size": str(comp_data.get("size", 1.0)),
                "units": "litre" # Example unit
            })
            # Store the SBML ID back into our compartment dict
            self.compartments[comp_kegg_id]['sbml_id'] = comp_sbml_id


        # --- Create Species ---
        listOfSpecies = ET.SubElement(model, "listOfSpecies")
        default_compartment_id = self.compartments.get('default', {}).get('sbml_id', list(self.compartment_sbml_ids)[0] if self.compartment_sbml_ids else 'default_comp_fallback') # Use first generated if default failed
        if not self.compartment_sbml_ids:
             logging.warning("No compartments defined or found, SBML species might be invalid without a compartment reference.")

        for species_kegg_id, species_data in self.species.items():
             if not species_data or 'name' not in species_data: # Skip incompletely parsed/failed entries
                 logging.warning(f"Skipping species {species_kegg_id} due to missing data.")
                 continue

             # Create a valid, unique SBML ID
             sbml_id = self._generate_unique_sbml_id(f"s_{species_kegg_id}", self.species_sbml_ids)
             species_data['sbml_id'] = sbml_id # Store for reaction referencing

             species = ET.SubElement(listOfSpecies, "species", {
                "id": sbml_id,
                "name": species_data.get("name", species_kegg_id).replace(';', ''), # Clean name a bit
                "compartment": default_compartment_id,
                "hasOnlySubstanceUnits": "false", # Assume concentration usually
                "boundaryCondition": "false", # Assume species change unless specified otherwise
                "constant": "false" # Assume variable concentration
             })

             # Add Annotations (KEGG ID, formula, charge)
             annotation = ET.SubElement(species, "annotation")
             rdf_rdf = ET.Element("rdf:RDF", {"xmlns:rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
                                             "xmlns:bqbiol": "http://biomodels.net/biology-qualifiers/"})
             desc = ET.SubElement(rdf_rdf, "rdf:Description", {"rdf:about": f"#{sbml_id}"})
             # Link to KEGG Compound entry
             is_version_of = ET.SubElement(desc, "bqbiol:is")
             bag = ET.SubElement(is_version_of, "rdf:Bag")
             ET.SubElement(bag, "rdf:li", {"rdf:resource": f"http://identifiers.org/kegg.compound/{species_kegg_id}"})
             annotation.append(rdf_rdf)

             # Add non-standard annotation for formula/charge if needed
             kegg_info = ET.SubElement(annotation, "kegg_info", {"xmlns": "http://www.example.com/kegg_info"}) # Custom namespace
             ET.SubElement(kegg_info, "kegg_id").text = species_kegg_id
             if species_data.get("formula"):
                 ET.SubElement(kegg_info, "formula").text = species_data["formula"]
             if species_data.get("charge") is not None:
                 ET.SubElement(kegg_info, "charge").text = str(species_data["charge"])
             if species_data.get("names"):
                 names_elem = ET.SubElement(kegg_info, "names")
                 for n in species_data["names"]:
                      ET.SubElement(names_elem, "name").text = n


        # --- Create Reactions ---
        listOfReactions = ET.SubElement(model, "listOfReactions")
        for reaction_kegg_id, reaction_data in self.reactions.items():
            if not reaction_data or 'substrates' not in reaction_data or 'products' not in reaction_data:
                 logging.warning(f"Skipping reaction {reaction_kegg_id} due to incomplete data (missing substrate/product list).")
                 continue
            # Skip reactions with no participants parsed
            if not reaction_data.get("substrates") and not reaction_data.get("products"):
                logging.warning(f"Skipping reaction {reaction_kegg_id} because no substrates or products were parsed.")
                continue

            # Create unique SBML ID
            rxn_sbml_id = self._generate_unique_sbml_id(f"r_{reaction_kegg_id}", self.reaction_sbml_ids)
            reaction_data['sbml_id'] = rxn_sbml_id

            reaction = ET.SubElement(listOfReactions, "reaction", {
                "id": rxn_sbml_id,
                "name": reaction_data.get("name", reaction_kegg_id).replace(';', ''),
                "reversible": "true" if reaction_data.get("reversible", False) else "false",
                "fast": "false" # Assume not fast equilibrium
            })

            # Add Annotations (KEGG ID, Enzymes)
            annotation = ET.SubElement(reaction, "annotation")
            rdf_rdf = ET.Element("rdf:RDF", {"xmlns:rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
                                            "xmlns:bqbiol": "http://biomodels.net/biology-qualifiers/"})
            desc = ET.SubElement(rdf_rdf, "rdf:Description", {"rdf:about": f"#{rxn_sbml_id}"})
            # Link to KEGG Reaction entry
            is_version_of = ET.SubElement(desc, "bqbiol:is")
            bag = ET.SubElement(is_version_of, "rdf:Bag")
            ET.SubElement(bag, "rdf:li", {"rdf:resource": f"http://identifiers.org/kegg.reaction/{reaction_kegg_id}"})
            # Link to Enzymes (if any)
            if reaction_data.get("enzymes"):
                 controlled_by = ET.SubElement(desc, "bqbiol:isVersionOf") # Relation enzyme -> reaction
                 bag_enz = ET.SubElement(controlled_by, "rdf:Bag")
                 for enzyme_id in reaction_data["enzymes"]:
                      # Try to guess if it's EC or KO for identifiers.org URL
                      if '.' in enzyme_id : # Assume EC number
                          ET.SubElement(bag_enz, "rdf:li", {"rdf:resource": f"http://identifiers.org/ec-code/{enzyme_id}"})
                      elif enzyme_id.startswith('K'): # Assume KO number
                           ET.SubElement(bag_enz, "rdf:li", {"rdf:resource": f"http://identifiers.org/kegg.orthology/{enzyme_id}"})
                      else:
                          logging.debug(f"Unrecognized enzyme format {enzyme_id} in reaction {reaction_kegg_id}, omitting identifier.org link.")
            annotation.append(rdf_rdf)

             # Add custom annotation for equation strings
            kegg_info = ET.SubElement(annotation, "kegg_info", {"xmlns": "http://www.example.com/kegg_info"}) # Custom namespace
            ET.SubElement(kegg_info, "kegg_id").text = reaction_kegg_id
            if reaction_data.get("definition_str"):
                 ET.SubElement(kegg_info, "definition").text = reaction_data["definition_str"]
            if reaction_data.get("equation_str"):
                 ET.SubElement(kegg_info, "equation").text = reaction_data["equation_str"]


            # --- Reactants (Substrates) ---
            if reaction_data.get("substrates"):
                listOfReactants = ET.SubElement(reaction, "listOfReactants")
                for comp_kegg_id, stoich in reaction_data["substrates"]:
                    if comp_kegg_id in self.species and 'sbml_id' in self.species[comp_kegg_id]:
                        species_sbml_id = self.species[comp_kegg_id]['sbml_id']
                        ET.SubElement(listOfReactants, "speciesReference", {
                            "species": species_sbml_id,
                            "stoichiometry": str(stoich),
                            "constant": "true" # Stoichiometry is constant
                        })
                    else:
                        logging.warning(f"Substrate {comp_kegg_id} in reaction {reaction_kegg_id} not found or has no SBML ID. Skipping in reaction list.")

            # --- Products ---
            if reaction_data.get("products"):
                listOfProducts = ET.SubElement(reaction, "listOfProducts")
                for comp_kegg_id, stoich in reaction_data["products"]:
                     if comp_kegg_id in self.species and 'sbml_id' in self.species[comp_kegg_id]:
                        species_sbml_id = self.species[comp_kegg_id]['sbml_id']
                        ET.SubElement(listOfProducts, "speciesReference", {
                            "species": species_sbml_id,
                            "stoichiometry": str(stoich),
                            "constant": "true"
                        })
                     else:
                         logging.warning(f"Product {comp_kegg_id} in reaction {reaction_kegg_id} not found or has no SBML ID. Skipping in reaction list.")


            # --- Kinetic Law (Generic Mass Action Placeholder) ---
            kineticLaw = ET.SubElement(reaction, "kineticLaw")
            math = ET.SubElement(kineticLaw, "math", {"xmlns": "http://www.w3.org/1998/Math/MathML"})

            # Parameters (k_forward, k_reverse)
            listOfParameters = ET.SubElement(kineticLaw, "listOfLocalParameters")
            param_k_fwd_id = self._generate_unique_sbml_id(f"kf_{rxn_sbml_id}", self.parameter_sbml_ids)
            ET.SubElement(listOfParameters, "localParameter", {"id": param_k_fwd_id, "value": "1.0", "units": "per_second"}) # Example units

            if reaction_data.get("reversible"):
                 param_k_rev_id = self._generate_unique_sbml_id(f"kr_{rxn_sbml_id}", self.parameter_sbml_ids)
                 ET.SubElement(listOfParameters, "localParameter", {"id": param_k_rev_id, "value": "0.1", "units": "per_second"})

            # Construct MathML for V = kf * [S1]^s1 * [S2]^s2 ... - kr * [P1]^p1 * [P2]^p2 ...
            # Or V = kf * [S1]^s1 * [S2]^s2 ... for irreversible

            forward_term = ET.Element("apply") # Represents kf * Product(reactants)
            ET.SubElement(forward_term, "times/>")
            ET.SubElement(forward_term, "ci").text = param_k_fwd_id
            for comp_kegg_id, stoich in reaction_data.get("substrates", []):
                if comp_kegg_id in self.species and 'sbml_id' in self.species[comp_kegg_id]:
                     species_sbml_id = self.species[comp_kegg_id]['sbml_id']
                     if stoich == 1:
                          ET.SubElement(forward_term, "ci").text = species_sbml_id
                     else:
                          power_apply = ET.SubElement(forward_term, "apply")
                          ET.SubElement(power_apply, "power/>")
                          ET.SubElement(power_apply, "ci").text = species_sbml_id
                          ET.SubElement(power_apply, "cn", {"type": "integer"}).text = str(stoich)

            if reaction_data.get("reversible"):
                reverse_term = ET.Element("apply") # Represents kr * Product(products)
                ET.SubElement(reverse_term, "times/>")
                ET.SubElement(reverse_term, "ci").text = param_k_rev_id
                for comp_kegg_id, stoich in reaction_data.get("products", []):
                    if comp_kegg_id in self.species and 'sbml_id' in self.species[comp_kegg_id]:
                         species_sbml_id = self.species[comp_kegg_id]['sbml_id']
                         if stoich == 1:
                             ET.SubElement(reverse_term, "ci").text = species_sbml_id
                         else:
                              power_apply = ET.SubElement(reverse_term, "apply")
                              ET.SubElement(power_apply, "power/>")
                              ET.SubElement(power_apply, "ci").text = species_sbml_id
                              ET.SubElement(power_apply, "cn", {"type": "integer"}).text = str(stoich)

                # Combine forward and reverse: V = forward - reverse
                main_apply = ET.SubElement(math, "apply")
                ET.SubElement(main_apply, "minus/>")
                main_apply.append(forward_term)
                main_apply.append(reverse_term)
            else:
                 # Just use forward term: V = forward
                 math.append(forward_term)


        logging.info("SBML model structure created.")
        return sbml


    def save_sbml(self, output_file: Optional[str] = None) -> Optional[str]:
        """
        Orchestrates parsing and saving the KEGG pathway to an SBML file.

        Args:
            output_file: Path to save the SBML file (default: {pathway_id}.xml).

        Returns:
            Path to the saved SBML file, or None if saving failed.
        """
        logging.info(f"Starting conversion process for KEGG pathway: {self.pathway_id}")

        # Reset internal state for potential reuse of the object
        self.species = {}
        self.reactions = {}
        self.species_sbml_ids = set()
        self.reaction_sbml_ids = set()
        self.parameter_sbml_ids = set()
        self.compartment_sbml_ids = set()
        self.compartments = {"default": {"name": "cytosol", "size": 1.0}} # Reset default compartment too


        # Parse the pathway - populates self.species and self.reactions
        success = self.parse_pathway()
        if not success:
             logging.error("Pathway parsing failed. Cannot generate SBML file.")
             return None
        if not self.species or not self.reactions:
             logging.error("Parsing seemed successful, but no species or reactions were collected. Cannot generate SBML file.")
             return None

        # Create the SBML model from the parsed data
        sbml_model_element = self.create_sbml_model()

        # Convert ElementTree to string with pretty formatting
        try:
            # Use minidom for pretty printing XML
            xml_str_bytes = ET.tostring(sbml_model_element, encoding='utf-8')
            dom = minidom.parseString(xml_str_bytes)
            pretty_xml = dom.toprettyxml(indent="  ", encoding='utf-8') # Specify encoding for output bytes
        except Exception as e:
            logging.error(f"Failed to format SBML XML for output: {e}")
             # Fallback to less pretty output if minidom fails
            pretty_xml = ET.tostring(sbml_model_element, encoding='utf-8', xml_declaration=True)

        # Determine output filename
        if output_file is None:
            output_file = f"{self.pathway_id}.xml"

        # Ensure directory exists if output_file includes a path
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                logging.info(f"Created output directory: {output_dir}")
            except OSError as e:
                 logging.error(f"Could not create output directory {output_dir}: {e}")
                 return None # Cannot save file

        # Save to file
        try:
            with open(output_file, 'wb') as f: # Write in binary mode due to encoding
                f.write(pretty_xml)
            logging.info(f"SBML file successfully saved as: {output_file}")
            return output_file
        except IOError as e:
            logging.error(f"Failed to write SBML file {output_file}: {e}")
            return None


def main():
    """Main function to run the converter from command line."""
    parser = argparse.ArgumentParser(
        description='Convert KEGG pathway data to SBML Level 3 Version 2 format.',
        epilog="Example: python kegg_to_sbml.py -p map00790 -o folate_biosynthesis.xml"
    )
    parser.add_argument(
        '-p', '--pathway',
        required=True, # Make pathway ID mandatory
        help='KEGG Pathway ID (e.g., map00790 for Folate Biosynthesis)'
    )
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output file name for the SBML document (default: {pathway_id}.xml)'
    )
    parser.add_argument(
        '-d', '--delay',
        type=float, default=0.1,
        help='Delay in seconds between KEGG API calls (default: 0.1)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging (DEBUG level)'
    )
    args = parser.parse_args()

    # Update logging level if verbose flag is set
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.debug("Verbose logging enabled.")

    try:
        converter = KEGGtoSBMLConverter(pathway_id=args.pathway, delay=args.delay)
        output_file = converter.save_sbml(output_file=args.output)

        if output_file:
            print(f"Conversion complete. SBML file saved as: {output_file}")
            # Simple validation: Check if file exists and is not empty
            if os.path.exists(output_file) and os.path.getsize(output_file) > 100: # Check size > 100 bytes as basic check
                print("Basic file validation successful (exists and not empty).")
            else:
                 print(f"Warning: Output file {output_file} might be empty or invalid.")
        else:
            print(f"Conversion failed for pathway {args.pathway}. See logs for details.")

    except ValueError as ve:
        print(f"Configuration Error: {ve}")
    except Exception as e:
         print(f"An unexpected error occurred: {e}")
         logging.exception("Unhandled exception during conversion process.")


if __name__ == "__main__":
    main()
