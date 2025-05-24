#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
import re
from typing import List, Tuple # Import List and Tuple for type hinting

class KGMLtoSBMLConverter:
    """Converts KGML (KEGG Markup Language) files to SBML format."""

    def __init__(self, kgml_filepath: str):
        """
        Initialize the converter.

        Args:
            kgml_filepath: Path to the input KGML file.
        """
        self.kgml_filepath = kgml_filepath
        self.kgml_tree = None
        self.kgml_root = None

        self.pathway_sbml_id = "default_pathway"
        self.pathway_name = "Default Pathway Name"
        self.pathway_organism = "" # e.g. "hsa"

        # Store parsed data:
        # self.species: {kegg_compound_id: {'sbml_id': '', 'name': '', 'kegg_ids': ['cpd:CXXXXX'], 'type': 'compound'}}
        #              {gene_entry_id: {'sbml_id': '', 'name': '', 'kegg_ids': ['hsa:X', 'ko:Y'], 'type': 'enzyme'}}
        self.species = {}
        
        # self.reactions: {kegg_reaction_id: {'sbml_id': '', 'name': '', 'reversible': False, 
        #                                     'substrates': [{'id': cpd_sbml_id, 'stoichiometry': 1}], 
        #                                     'products': [{'id': cpd_sbml_id, 'stoichiometry': 1}],
        #                                     'modifiers': [{'id': enzyme_sbml_id}] }}
        self.reactions = {}

        # For managing unique SBML IDs
        self.sbml_id_set = set() # Stores all generated SBML IDs to ensure uniqueness

    def _get_unique_sbml_id(self, base_id: str, prefix: str = "") -> str:
        """Generates a unique SBML ID."""
        # Sanitize base_id: remove invalid characters, ensure it starts with a letter or underscore
        processed_base_id = re.sub(r'[^A-Za-z0-9_]', '_', base_id)
        if not re.match(r'^[A-Za-z_]', processed_base_id):
            processed_base_id = f"_{processed_base_id}"
        
        full_base_id = f"{prefix}{processed_base_id}"
        
        sbml_id = full_base_id
        n = 1
        while sbml_id in self.sbml_id_set:
            sbml_id = f"{full_base_id}_{n}"
            n += 1
        self.sbml_id_set.add(sbml_id)
        return sbml_id

    def _sanitize_kegg_id_for_uri(self, kegg_id_str: str) -> str:
        """Extracts the core ID from a KEGG ID string like 'cpd:C00001' -> 'C00001'."""
        if ':' in kegg_id_str:
            return kegg_id_str.split(':')[-1]
        return kegg_id_str
        
    def _extract_kegg_ids_from_name_attr(self, name_attr: str) -> list:
        """Extracts multiple KEGG IDs from an entry's name attribute (e.g., "hsa:123 hsa:456")"""
        return name_attr.split()

    def _parse_kgml(self):
        """Parses the KGML file and populates internal data structures."""
        try:
            self.kgml_tree = ET.parse(self.kgml_filepath)
            self.kgml_root = self.kgml_tree.getroot()
        except ET.ParseError as e:
            raise ValueError(f"Error parsing KGML file: {e}")
        except FileNotFoundError:
            raise FileNotFoundError(f"KGML file not found: {self.kgml_filepath}")

        # 1. Parse Pathway Information
        pathway_kegg_id = self.kgml_root.get('name') # e.g. path:hsa00010
        self.pathway_sbml_id = self._get_unique_sbml_id(self._sanitize_kegg_id_for_uri(pathway_kegg_id), "model_")
        self.pathway_name = self.kgml_root.get('title', 'Untitled Pathway')
        self.pathway_organism = self.kgml_root.get('org', '')

        # Temporary mapping for gene entries to their associated reaction KEGG IDs
        gene_entry_to_reactions = {} # {gene_entry_kgml_id: [rn_id1, rn_id2]}

        # 2. First Pass: Parse Entries (compounds, genes/enzymes)
        kgml_entries = {} # {kgml_entry_id: entry_data}
        for entry_elem in self.kgml_root.findall('entry'):
            entry_id = entry_elem.get('id')
            entry_type = entry_elem.get('type')
            # KEGG name attribute (e.g., "cpd:C00001", "hsa:123 ko:K456")
            entry_kegg_name_attr = entry_elem.get('name') 
            
            graphics_elem = entry_elem.find('graphics')
            # Name from graphics, often more human-readable or primary name
            graphics_name_str = graphics_elem.get('name', '') if graphics_elem is not None else ''
            # Clean up "..." if present
            graphics_name_str = graphics_name_str.replace('...', '').strip()

            kgml_entries[entry_id] = {
                'kegg_name_attr': entry_kegg_name_attr,
                'type': entry_type,
                'graphics_name': graphics_name_str
            }

            if entry_type == 'compound':
                compound_kegg_id = self._sanitize_kegg_id_for_uri(entry_kegg_name_attr) # e.g. C00001
                if compound_kegg_id not in self.species:
                    sbml_id = self._get_unique_sbml_id(compound_kegg_id, "s_")
                    self.species[compound_kegg_id] = {
                        'sbml_id': sbml_id,
                        'name': graphics_name_str if graphics_name_str else compound_kegg_id,
                        'kegg_ids_full': self._extract_kegg_ids_from_name_attr(entry_kegg_name_attr), # ['cpd:CXXXX']
                        'type': 'compound'
                    }
                # Store SBML ID back into kgml_entries for cross-referencing if needed
                kgml_entries[entry_id]['sbml_id'] = self.species[compound_kegg_id]['sbml_id']
            
            elif entry_type == 'gene' or entry_type == 'ortholog':
                # For genes/orthologs, we use the KGML entry ID as the key for self.species dict,
                # as one entry might represent multiple gene products acting together or as isoenzymes.
                if entry_id not in self.species:
                    base_id_for_sbml = graphics_name_str.split(',')[0].strip() # Use first name from graphics if available
                    if not base_id_for_sbml: # Fallback if graphics name is empty
                        base_id_for_sbml = self._sanitize_kegg_id_for_uri(entry_kegg_name_attr.split()[0]) # Use first KEGG ID
                    
                    sbml_id = self._get_unique_sbml_id(f"enzyme_{entry_id}_{base_id_for_sbml}", "s_")
                    
                    display_name = graphics_name_str if graphics_name_str else entry_kegg_name_attr

                    self.species[entry_id] = { # Using entry_id as key for enzymes for now
                        'sbml_id': sbml_id,
                        'name': display_name,
                        'kegg_ids_full': self._extract_kegg_ids_from_name_attr(entry_kegg_name_attr), # ['hsa:X', 'ko:Y']
                        'type': 'enzyme' # Generic type for SBML species
                    }
                kgml_entries[entry_id]['sbml_id'] = self.species[entry_id]['sbml_id']

                # Store reaction associations for this gene/ortholog entry
                reaction_attr = entry_elem.get('reaction')
                if reaction_attr:
                    gene_entry_to_reactions[entry_id] = [self._sanitize_kegg_id_for_uri(r) for r in reaction_attr.split()]


        # 3. Second Pass: Parse Reactions from KGML <reaction> elements
        for reaction_elem in self.kgml_root.findall('reaction'):
            reaction_kegg_id_full = reaction_elem.get('name') # e.g. rn:R00001
            reaction_kegg_id = self._sanitize_kegg_id_for_uri(reaction_kegg_id_full) # R00001
            
            if reaction_kegg_id in self.reactions: # Already processed (e.g. from another definition)
                continue

            sbml_reaction_id = self._get_unique_sbml_id(reaction_kegg_id, "r_")
            reaction_name = reaction_kegg_id # Default name
            reaction_reversible = reaction_elem.get('type') == 'reversible'
            
            current_reaction_data = {
                'sbml_id': sbml_reaction_id,
                'name': reaction_name,
                'reversible': reaction_reversible,
                'substrates': [],
                'products': [],
                'modifiers': [],
                'kegg_id_full': reaction_kegg_id_full
            }

            # Parse substrates
            for substrate_elem in reaction_elem.findall('substrate'):
                # The 'name' attribute of substrate/product in KGML is the KEGG compound ID
                compound_kegg_id_full = substrate_elem.get('name') # e.g., cpd:C00001
                compound_kegg_id = self._sanitize_kegg_id_for_uri(compound_kegg_id_full)
                
                if compound_kegg_id in self.species and self.species[compound_kegg_id]['type'] == 'compound':
                    substrate_sbml_id = self.species[compound_kegg_id]['sbml_id']
                    current_reaction_data['substrates'].append({'id': substrate_sbml_id, 'stoichiometry': 1})
                else:
                    # This compound wasn't defined as an <entry type="compound">. Create it.
                    if compound_kegg_id not in self.species:
                        s_id = self._get_unique_sbml_id(compound_kegg_id, "s_")
                        self.species[compound_kegg_id] = {
                            'sbml_id': s_id, 'name': compound_kegg_id,
                            'kegg_ids_full': [compound_kegg_id_full], 'type': 'compound'
                        }
                        current_reaction_data['substrates'].append({'id': s_id, 'stoichiometry': 1})
                    elif self.species[compound_kegg_id]['type'] == 'compound': # exists but check if missed linking somehow
                         current_reaction_data['substrates'].append({'id': self.species[compound_kegg_id]['sbml_id'], 'stoichiometry': 1})


            # Parse products
            for product_elem in reaction_elem.findall('product'):
                compound_kegg_id_full = product_elem.get('name') # e.g., cpd:C00002
                compound_kegg_id = self._sanitize_kegg_id_for_uri(compound_kegg_id_full)

                if compound_kegg_id in self.species and self.species[compound_kegg_id]['type'] == 'compound':
                    product_sbml_id = self.species[compound_kegg_id]['sbml_id']
                    current_reaction_data['products'].append({'id': product_sbml_id, 'stoichiometry': 1})
                else:
                    if compound_kegg_id not in self.species:
                        s_id = self._get_unique_sbml_id(compound_kegg_id, "s_")
                        self.species[compound_kegg_id] = {
                            'sbml_id': s_id, 'name': compound_kegg_id,
                            'kegg_ids_full': [compound_kegg_id_full], 'type': 'compound'
                        }
                        current_reaction_data['products'].append({'id': s_id, 'stoichiometry': 1})
                    elif self.species[compound_kegg_id]['type'] == 'compound':
                         current_reaction_data['products'].append({'id': self.species[compound_kegg_id]['sbml_id'], 'stoichiometry': 1})

            self.reactions[reaction_kegg_id] = current_reaction_data

        # 4. Link Enzymes (from gene entries) to Reactions as Modifiers
        #    And potentially update reaction names
        enzyme_names_for_reactions = {} # {reaction_kegg_id: [enzyme_name1, ...]}

        for gene_entry_id, associated_rn_ids in gene_entry_to_reactions.items():
            if gene_entry_id in self.species and self.species[gene_entry_id]['type'] == 'enzyme':
                enzyme_sbml_id = self.species[gene_entry_id]['sbml_id']
                enzyme_display_name = self.species[gene_entry_id]['name']
                for rn_id in associated_rn_ids:
                    if rn_id in self.reactions:
                        self.reactions[rn_id]['modifiers'].append({'id': enzyme_sbml_id})
                        if rn_id not in enzyme_names_for_reactions:
                            enzyme_names_for_reactions[rn_id] = []
                        if enzyme_display_name:
                             enzyme_names_for_reactions[rn_id].append(enzyme_display_name)
        
        # Update reaction names with enzyme names if available
        for rn_id, enzyme_names in enzyme_names_for_reactions.items():
            if rn_id in self.reactions and enzyme_names:
                # Only use unique enzyme names for the reaction name
                unique_enz_names = sorted(list(set(name.split(',')[0].strip() for name in enzyme_names if name))) # Take first part of name
                if unique_enz_names:
                    self.reactions[rn_id]['name'] = f"{rn_id}: {', '.join(unique_enz_names)}"


    def _create_sbml_annotation(self, parent_element: ET.Element, sbml_element_id: str, kegg_id_tuples: List[Tuple[str, str]]):
        """
        Creates RDF annotation for an SBML element.
        kegg_id_tuples: list of (kegg_type_prefix, kegg_id_full) e.g., [('compound', 'cpd:C00001')] or [('genes', 'hsa:123'), ('orthology', 'ko:K456')]
        """
        if not kegg_id_tuples:
            return

        annotation = ET.SubElement(parent_element, "annotation")
        rdf_rdf = ET.SubElement(annotation, "rdf:RDF", {
            "xmlns:rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "xmlns:bqbiol": "http://biomodels.net/biology-qualifiers/"
        })
        description = ET.SubElement(rdf_rdf, "rdf:Description", {"rdf:about": f"#{sbml_element_id}"})
        
        # Use bqbiol:is for single primary ID, bqbiol:isVersionOf or bqbiol:hasPart for multiple related IDs
        # For simplicity here, we'll use bqbiol:is for all, linking to each provided KEGG ID.
        
        # Determine the correct BioModels qualifier
        # Assuming the first ID in kegg_id_tuples defines the "primary" type for bqbiol:is or bqbiol:isVersionOf
        # If multiple distinct KEGG IDs are present (e.g. for a gene entry `name="hsa:X ko:Y"`), `isVersionOf` might be better.
        # If it's just one ID (e.g. a compound or a reaction), `is` is fine.
        qualifier_tag = "bqbiol:is"
        # A more robust check for multiple distinct types or a complex entity
        is_complex_entity = False
        if len(kegg_id_tuples) > 1:
            first_type = kegg_id_tuples[0][0]
            # Check if there are other types or if any ID represents a collection (like 'genes' for multiple gene products)
            if any(t[0] != first_type for t in kegg_id_tuples) or \
               any(t[0] in ["genes", "orthology"] and " " in t[1] for t in kegg_id_tuples): # Simplistic check for multiple items in one ID string
                is_complex_entity = True
        
        if is_complex_entity:
            qualifier_tag = "bqbiol:isDescribedBy" # Using for more complex cases like gene groups

        
        bq_element = ET.SubElement(description, qualifier_tag)
        bag = ET.SubElement(bq_element, "rdf:Bag")

        for kegg_db_prefix, full_kegg_id in kegg_id_tuples:
            # kegg_db_prefix like "compound", "reaction", "genes", "orthology"
            # full_kegg_id like "cpd:CXXXX", "rn:RXXXX", "hsa:123", "ko:K123"
            
            # Split if full_kegg_id itself contains multiple space-separated IDs (e.g. "hsa:123 hsa:456")
            # This can happen if _extract_kegg_ids_from_name_attr wasn't fully utilized upstream for this specific tuple generation
            individual_ids = full_kegg_id.split()

            for single_id in individual_ids:
                core_id = self._sanitize_kegg_id_for_uri(single_id)
                uri = ""

                if kegg_db_prefix == "compound":
                    uri = f"http://identifiers.org/kegg.compound/{core_id}"
                elif kegg_db_prefix == "reaction":
                    uri = f"http://identifiers.org/kegg.reaction/{core_id}"
                elif kegg_db_prefix == "genes": # e.g., hsa:123, eco:b0001
                    # The full_kegg_id here should be like "hsa:123"
                    # identifiers.org format is kegg.genes/<org_code>:<gene_id> or kegg.genes/<entry>
                    # However, KEGG API links often use kegg.genes/ko:Knumber or kegg.genes/hsa:123
                    # Let's assume single_id is like "hsa:123" or "eco:b0001"
                    uri = f"http://identifiers.org/kegg.genes/{single_id}"
                elif kegg_db_prefix == "orthology": # e.g., ko:K12345
                    uri = f"http://identifiers.org/kegg.orthology/{core_id}" # core_id would be K12345
                elif kegg_db_prefix.startswith("path"): # e.g. path:hsa, path:ko. full_kegg_id = "path:hsa00010"
                    # The core_id will be e.g., "hsa00010" or "ko00010"
                    uri = f"http://identifiers.org/kegg.pathway/{core_id}"
                else: # Fallback for unknown or less common types
                    # Try to construct a generic one if possible, or skip
                    if ":" in single_id:
                        db, actual_id = single_id.split(":", 1)
                        uri = f"http://identifiers.org/kegg.{db}/{actual_id}" # Best guess
                    else:
                        # Cannot form a good URI, skip this particular ID
                        print(f"Warning: Could not form Identifiers.org URI for prefix '{kegg_db_prefix}' and ID '{single_id}' for SBML element '{sbml_element_id}'. Skipping annotation for this ID.")
                        continue
                
                if uri:
                    ET.SubElement(bag, "rdf:li", {"rdf:resource": uri})


    def create_sbml_model(self) -> ET.Element:
        """Creates the SBML ET.Element tree from parsed KGML data."""
        sbml_ns = "http://www.sbml.org/sbml/level3/version2/core" # L3V2 is common
        sbml = ET.Element("sbml", {"xmlns": sbml_ns, "level": "3", "version": "2"})
        model = ET.SubElement(sbml, "model", {"id": self.pathway_sbml_id, "name": self.pathway_name})

        # Notes
        notes_str = f"Model generated from KEGG KGML file '{self.kgml_filepath}' for pathway {self.pathway_name} ({self._sanitize_kegg_id_for_uri(self.kgml_root.get('name'))}). Conversion performed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}."
        notes_elem = ET.SubElement(model, "notes")
        html_notes = ET.SubElement(notes_elem, "html", {"xmlns": "http://www.w3.org/1999/xhtml"})
        p_elem = ET.SubElement(html_notes, "p")
        p_elem.text = notes_str
        
        # Pathway annotation for the model itself
        model_kegg_ids_tuples = []
        pathway_kegg_full_id = self.kgml_root.get('name') # e.g. path:hsa00010
        if pathway_kegg_full_id:
             # The prefix should represent the database, for pathway it's 'pathway' but it also contains organism info
             # For identifiers.org, it's kegg.pathway/map00010 or kegg.pathway/hsa00010
             # The 'kegg_db_prefix' is used to construct identifiers.org/kegg.{prefix}/
             # So we need to pass the actual ID like "hsa00010" and a prefix like "pathway"
             # The _sanitize_kegg_id_for_uri will strip 'path:' part.
            sanitized_path_id = self._sanitize_kegg_id_for_uri(pathway_kegg_full_id) # e.g. hsa00010
            model_kegg_ids_tuples.append(("pathway", sanitized_path_id)) # Use "pathway" as prefix, pass full "hsa00010"
        self._create_sbml_annotation(model, self.pathway_sbml_id, model_kegg_ids_tuples)


        # Compartments
        listOfCompartments = ET.SubElement(model, "listOfCompartments")
        default_comp_id = "default_compartment"
        ET.SubElement(listOfCompartments, "compartment", {
            "id": default_comp_id, "name": "Default Compartment", 
            "constant": "true", "spatialDimensions": "3", "size": "1"
        })

        # Species
        listOfSpecies = ET.SubElement(model, "listOfSpecies")
        # self.species keys can be compound_kegg_ids OR gene_entry_ids
        for _, species_data in self.species.items():
            sbml_id = species_data['sbml_id']
            name = species_data['name']
            kegg_ids_full_list = species_data['kegg_ids_full'] # list of 'cpd:CXX', 'hsa:X', 'ko:KXX'
            species_type = species_data['type'] # 'compound' or 'enzyme'

            species_elem = ET.SubElement(listOfSpecies, "species", {
                "id": sbml_id, "name": name, "compartment": default_comp_id,
                "hasOnlySubstanceUnits": "false", "boundaryCondition": "false", "constant": "false"
            })
            
            # Annotations for species
            kegg_id_tuples_for_annotation = []
            for kid_full in kegg_ids_full_list: # kid_full is e.g. "cpd:C00001", "hsa:123", "ko:K456"
                db_prefix = "unknown"
                if kid_full.startswith("cpd:"): db_prefix = "compound"
                elif kid_full.startswith("gl:"): db_prefix = "glycan" # kegg.glycan
                elif kid_full.startswith("ko:"): db_prefix = "orthology" # kegg.orthology
                elif re.match(r"^[a-z]{3,4}:", kid_full): # e.g., hsa:, eco: (KEGG Gene ID)
                    db_prefix = "genes" # It refers to kegg.genes/<org_id_pair> e.g. kegg.genes/hsa:123
                
                if db_prefix != "unknown":
                    kegg_id_tuples_for_annotation.append((db_prefix, kid_full)) # Pass the full ID "hsa:123"
                else:
                    print(f"Warning: Unknown KEGG ID prefix for '{kid_full}' in species '{name}'. Skipping annotation for this ID.")

            self._create_sbml_annotation(species_elem, sbml_id, kegg_id_tuples_for_annotation)


        # Reactions
        listOfReactions = ET.SubElement(model, "listOfReactions")
        for _, reaction_data in self.reactions.items():
            sbml_id = reaction_data['sbml_id']
            name = reaction_data['name']
            kegg_reaction_id_full = reaction_data['kegg_id_full'] # e.g. rn:R00001

            reaction_elem = ET.SubElement(listOfReactions, "reaction", {
                "id": sbml_id, "name": name, 
                "reversible": "true" if reaction_data['reversible'] else "false",
                "fast": "false" # Default
            })
            
            # The full ID like "rn:R12345" is passed, _sanitize_kegg_id_for_uri will extract "R12345" in _create_sbml_annotation
            self._create_sbml_annotation(reaction_elem, sbml_id, [('reaction', kegg_reaction_id_full)])


            # Reactants
            if reaction_data['substrates']:
                listOfReactants = ET.SubElement(reaction_elem, "listOfReactants")
                for s_info in reaction_data['substrates']:
                    ET.SubElement(listOfReactants, "speciesReference", {
                        "species": s_info['id'], "stoichiometry": str(s_info['stoichiometry']), "constant": "true"
                    })
            
            # Products
            if reaction_data['products']:
                listOfProducts = ET.SubElement(reaction_elem, "listOfProducts")
                for p_info in reaction_data['products']:
                    ET.SubElement(listOfProducts, "speciesReference", {
                        "species": p_info['id'], "stoichiometry": str(p_info['stoichiometry']), "constant": "true"
                    })

            # Modifiers (Enzymes)
            if reaction_data['modifiers']:
                listOfModifiers = ET.SubElement(reaction_elem, "listOfModifiers")
                for m_info in reaction_data['modifiers']:
                    ET.SubElement(listOfModifiers, "modifierSpeciesReference", {"species": m_info['id']})
            
            # Kinetic Law (default mass action)
            kineticLaw = ET.SubElement(reaction_elem, "kineticLaw")
            # MathML
            math = ET.SubElement(kineticLaw, "math", {"xmlns": "http://www.w3.org/1998/Math/MathML"})
            
            # Create kinetic law string: k_f * S1 * S2 ... or k_f * S1 * S2 - k_r * P1 * P2 ...
            # Parameters for kinetic law
            listOfLocalParameters = ET.SubElement(kineticLaw, "listOfLocalParameters") # L3V2 uses listOfLocalParameters
            
            k_forward_id = self._get_unique_sbml_id(f"kf_{sbml_id}")
            ET.SubElement(listOfLocalParameters, "localParameter", {"id": k_forward_id, "value": "0.1", "units": "dimensionless"})
            
            term_parts = [f'<ci>{k_forward_id}</ci>']
            for reactant in reaction_data['substrates']:
                term_parts.append(f'<ci>{reactant["id"]}</ci>')
            
            forward_term_str = ""
            if len(term_parts) > 1: # Has kf and at least one reactant
                 forward_term_str = "<apply><times/>" + "".join(term_parts) + "</apply>"
            elif len(term_parts) == 1 and not reaction_data['substrates']: # only kf, zero-order forward
                 forward_term_str = term_parts[0]

            full_math_str = forward_term_str

            if reaction_data['reversible']:
                k_reverse_id = self._get_unique_sbml_id(f"kr_{sbml_id}")
                ET.SubElement(listOfLocalParameters, "localParameter", {"id": k_reverse_id, "value": "0.01", "units": "dimensionless"})
                
                rev_term_parts = [f'<ci>{k_reverse_id}</ci>']
                for product in reaction_data['products']:
                    rev_term_parts.append(f'<ci>{product["id"]}</ci>')

                reverse_term_str = ""
                if len(rev_term_parts) > 1 : # Has kr and at least one product
                    reverse_term_str = "<apply><times/>" + "".join(rev_term_parts) + "</apply>"
                elif len(rev_term_parts) == 1 and not reaction_data['products']: # only kr, zero-order reverse
                    reverse_term_str = rev_term_parts[0]

                if forward_term_str and reverse_term_str:
                    full_math_str = f"<apply><minus/>{forward_term_str}{reverse_term_str}</apply>"
                elif forward_term_str: # Only forward part if reverse is empty
                    full_math_str = forward_term_str
                elif reverse_term_str: # Only reverse part if forward is empty (unusual but possible if no reactants)
                     full_math_str = f"<apply><minus/><cn type='real'>0</cn>{reverse_term_str}</apply>" # 0 - k_r * P


            # For reactions with no reactants and no products (e.g. boundary species exchange)
            # or just one direction having no species
            if not full_math_str:
                 full_math_str = f"<cn type='real'>0.0</cn>" # Default to zero rate if expression is empty


            # Parse the string and add to math element - this is a hack for simplicity
            # A proper MathML builder should be used for complex cases.
            try:
                math_expr_elem = ET.fromstring(full_math_str)
                math.append(math_expr_elem)
            except ET.ParseError: # Fallback if string is not valid XML
                cn_elem = ET.SubElement(math, "cn", type="real") # Default to 0 if math string parsing fails
                cn_elem.text = "0.0"
                print(f"Warning: Could not parse MathML string for reaction {sbml_id}: {full_math_str}. Defaulting to rate 0.")

        return sbml


    def convert_and_save(self, output_filepath: str):
        """
        Parses the KGML, converts to SBML, and saves to the output file.

        Args:
            output_filepath: Path to save the generated SBML file.
        """
        print(f"Parsing KGML file: {self.kgml_filepath}")
        self._parse_kgml()
        
        print("Generating SBML model...")
        sbml_model_element = self.create_sbml_model()
        
        # Pretty print XML
        xml_str = ET.tostring(sbml_model_element, encoding='utf-8', method='xml')
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ")
        
        with open(output_filepath, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)
        print(f"SBML file saved as: {output_filepath}")


def main():
    parser = argparse.ArgumentParser(description='Convert KEGG KGML file to SBML format.')
    parser.add_argument('-i', '--input', required=True, help='Input KGML file path.')
    parser.add_argument('-o', '--output', required=True, help='Output SBML file path (e.g., model.xml).')
    args = parser.parse_args()

    try:
        converter = KGMLtoSBMLConverter(args.input)
        converter.convert_and_save(args.output)
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
