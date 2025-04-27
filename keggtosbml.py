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

class KEGGtoSBMLConverter:
    """Convert KEGG pathway data to SBML format."""
    
    def __init__(self, pathway_id: str = "map00790"):
        """
        Initialize the converter with a KEGG pathway ID.
        
        Args:
            pathway_id: KEGG pathway ID (default: map00790 - Folate Biosynthesis)
        """
        self.pathway_id = pathway_id
        self.base_url = "http://rest.kegg.jp"
        self.species = {}  # Dictionary to store species/compounds
        self.reactions = {}  # Dictionary to store reactions
        self.pathway_name = ""
        self.compartments = {"default": "cytosol"}
        
        # For managing SBML elements
        self.species_ids = set()
        self.reaction_ids = set()
        
    def fetch_pathway_data(self) -> str:
        """
        Fetch pathway data from KEGG API.
        
        Returns:
            Raw KEGG pathway data as string
        """
        url = f"{self.base_url}/get/{self.pathway_id}"
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception(f"Failed to fetch pathway data: {response.status_code}")
        return response.text
    
    def fetch_pathway_name(self) -> str:
        """
        Fetch the name of the pathway.
        
        Returns:
            Pathway name
        """
        url = f"{self.base_url}/find/pathway/{self.pathway_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return f"Pathway {self.pathway_id}"
        
        # Extract name from the response
        lines = response.text.strip().split('\n')
        if lines:
            parts = lines[0].split('\t')
            if len(parts) > 1:
                return parts[1].split(' - ')[0].strip()
        return f"Pathway {self.pathway_id}"
    
    def parse_compound_data(self, compound_id: str) -> Dict:
        """
        Fetch and parse compound data.
        
        Args:
            compound_id: KEGG compound ID (e.g., C00001)
            
        Returns:
            Dictionary with compound information
        """
        url = f"{self.base_url}/get/{compound_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return {"name": compound_id, "formula": "", "charge": 0}
        
        data = response.text
        compound_info = {"id": compound_id, "name": compound_id, "formula": "", "charge": 0}
        
        # Extract name
        name_match = re.search(r"NAME\s+(.*?)(?:\n\w+|\Z)", data, re.DOTALL)
        if name_match:
            compound_info["name"] = name_match.group(1).strip().split('\n')[0].strip()
        
        # Extract formula
        formula_match = re.search(r"FORMULA\s+(.*?)(?:\n\w+|\Z)", data, re.DOTALL)
        if formula_match:
            compound_info["formula"] = formula_match.group(1).strip()
        
        # Extract charge if available
        charge_match = re.search(r"CHARGE\s+([-+]?\d+)", data)
        if charge_match:
            compound_info["charge"] = int(charge_match.group(1))
            
        return compound_info
    
    def parse_reaction_data(self, reaction_id: str) -> Dict:
        """
        Fetch and parse reaction data.
        
        Args:
            reaction_id: KEGG reaction ID (e.g., R00001)
            
        Returns:
            Dictionary with reaction information
        """
        url = f"{self.base_url}/get/{reaction_id}"
        response = requests.get(url)
        if response.status_code != 200:
            return {"name": reaction_id, "equation": ""}
        
        data = response.text
        reaction_info = {"id": reaction_id, "name": reaction_id, "equation": "", 
                         "substrates": [], "products": [], "reversible": False}
        
        # Extract name
        name_match = re.search(r"NAME\s+(.*?)(?:\n\w+|\Z)", data, re.DOTALL)
        if name_match:
            reaction_info["name"] = name_match.group(1).strip().split('\n')[0].strip()
        
        # Extract equation
        equation_match = re.search(r"EQUATION\s+(.*?)(?:\n\w+|\Z)", data, re.DOTALL)
        if equation_match:
            equation = equation_match.group(1).strip()
            reaction_info["equation"] = equation
            
            # Determine reversibility
            reaction_info["reversible"] = "<=>" in equation
            
            # Parse substrates and products from equation
            if "<=>" in equation:
                substrates_str, products_str = equation.split("<=>")
            elif "=>" in equation:
                substrates_str, products_str = equation.split("=>")
            else:
                substrates_str, products_str = "", equation
                
            # Process substrates
            for substrate in substrates_str.split("+"):
                substrate = substrate.strip()
                if not substrate:
                    continue
                    
                # Check for stoichiometry
                stoich_match = re.match(r"(\d+)\s+(.+)", substrate)
                if stoich_match:
                    stoich = int(stoich_match.group(1))
                    comp_id = stoich_match.group(2).strip()
                else:
                    stoich = 1
                    comp_id = substrate.strip()
                    
                # Extract compound ID from format like "C00001" or "n C00001"
                compound_match = re.search(r"(C\d+)", comp_id)
                if compound_match:
                    comp_id = compound_match.group(1)
                    reaction_info["substrates"].append((comp_id, stoich))
                    
            # Process products
            for product in products_str.split("+"):
                product = product.strip()
                if not product:
                    continue
                    
                # Check for stoichiometry
                stoich_match = re.match(r"(\d+)\s+(.+)", product)
                if stoich_match:
                    stoich = int(stoich_match.group(1))
                    comp_id = stoich_match.group(2).strip()
                else:
                    stoich = 1
                    comp_id = product.strip()
                    
                # Extract compound ID
                compound_match = re.search(r"(C\d+)", comp_id)
                if compound_match:
                    comp_id = compound_match.group(1)
                    reaction_info["products"].append((comp_id, stoich))
                    
        return reaction_info
    
    def parse_pathway(self) -> None:
        """Parse the KEGG pathway data."""
        self.pathway_name = self.fetch_pathway_name()
        pathway_data = self.fetch_pathway_data()
        
        # Process compounds and reactions
        compound_pattern = r"COMPOUND\s+(.*?)(?:\n\w+|\Z)"
        compound_match = re.search(compound_pattern, pathway_data, re.DOTALL)
        if compound_match:
            compounds_text = compound_match.group(1).strip()
            for line in compounds_text.split('\n'):
                if line.strip():
                    parts = line.strip().split()
                    if parts and parts[0].startswith('C'):
                        compound_id = parts[0]
                        if compound_id not in self.species:
                            self.species[compound_id] = self.parse_compound_data(compound_id)
        
        # Parse reactions from REACTION section
        reaction_pattern = r"REACTION\s+(.*?)(?:\n\w+|\Z)"
        reaction_match = re.search(reaction_pattern, pathway_data, re.DOTALL)
        if reaction_match:
            reactions_text = reaction_match.group(1).strip()
            for line in reactions_text.split('\n'):
                if line.strip():
                    parts = line.strip().split()
                    if parts and parts[0].startswith('R'):
                        reaction_id = parts[0]
                        if reaction_id not in self.reactions:
                            reaction_data = self.parse_reaction_data(reaction_id)
                            self.reactions[reaction_id] = reaction_data
                            
                            # Add any new species found in reactions
                            for comp_id, _ in reaction_data["substrates"] + reaction_data["products"]:
                                if comp_id not in self.species:
                                    self.species[comp_id] = self.parse_compound_data(comp_id)
    
    def create_sbml_model(self) -> ET.Element:
        """
        Create SBML model from parsed KEGG data.
        
        Returns:
            XML Element containing the SBML model
        """
        # Create SBML structure
        sbml = ET.Element("sbml", {
            "xmlns": "http://www.sbml.org/sbml/level3/version1/core",
            "level": "3",
            "version": "1"
        })
        
        # Create model
        model = ET.SubElement(sbml, "model", {"id": f"KEGG_{self.pathway_id}"})
        ET.SubElement(model, "name").text = self.pathway_name
        
        # Add notes about the model creation
        notes = ET.SubElement(model, "notes")
        html_notes = ET.SubElement(notes, "html", {"xmlns": "http://www.w3.org/1999/xhtml"})
        body = ET.SubElement(html_notes, "body")
        ET.SubElement(body, "p").text = f"This model was automatically generated from KEGG pathway {self.pathway_id} ({self.pathway_name}) on {datetime.now().strftime('%Y-%m-%d')}."
        
        # Create compartments list
        listOfCompartments = ET.SubElement(model, "listOfCompartments")
        for comp_id, comp_name in self.compartments.items():
            compartment = ET.SubElement(listOfCompartments, "compartment", {
                "id": comp_id,
                "name": comp_name,
                "constant": "true",
                "spatialDimensions": "3",
                "size": "1"
            })
        
        # Create species list
        listOfSpecies = ET.SubElement(model, "listOfSpecies")
        for species_id, species_data in self.species.items():
            # Create a valid SBML ID
            sbml_id = f"s_{species_id}"
            n = 1
            while sbml_id in self.species_ids:
                sbml_id = f"s_{species_id}_{n}"
                n += 1
                
            self.species_ids.add(sbml_id)
            species_data['sbml_id'] = sbml_id
            
            species = ET.SubElement(listOfSpecies, "species", {
                "id": sbml_id,
                "name": species_data.get("name", species_id),
                "compartment": "default",
                "hasOnlySubstanceUnits": "false",
                "boundaryCondition": "false",
                "constant": "false"
            })
            
            # Add formula as annotation if available
            if species_data.get("formula"):
                annotation = ET.SubElement(species, "annotation")
                ET.SubElement(annotation, "formula").text = species_data["formula"]
            
        # Create reactions list
        listOfReactions = ET.SubElement(model, "listOfReactions")
        for reaction_id, reaction_data in self.reactions.items():
            # Create a valid SBML ID
            sbml_id = f"r_{reaction_id}"
            n = 1
            while sbml_id in self.reaction_ids:
                sbml_id = f"r_{reaction_id}_{n}"
                n += 1
                
            self.reaction_ids.add(sbml_id)
            
            reaction = ET.SubElement(listOfReactions, "reaction", {
                "id": sbml_id,
                "name": reaction_data.get("name", reaction_id),
                "reversible": "true" if reaction_data.get("reversible", False) else "false",
                "fast": "false"
            })
            
            # Add reactants (substrates)
            if reaction_data.get("substrates"):
                listOfReactants = ET.SubElement(reaction, "listOfReactants")
                for comp_id, stoich in reaction_data["substrates"]:
                    if comp_id in self.species:
                        species_ref = ET.SubElement(listOfReactants, "speciesReference", {
                            "species": self.species[comp_id]["sbml_id"],
                            "stoichiometry": str(stoich),
                            "constant": "true"
                        })
            
            # Add products
            if reaction_data.get("products"):
                listOfProducts = ET.SubElement(reaction, "listOfProducts")
                for comp_id, stoich in reaction_data["products"]:
                    if comp_id in self.species:
                        species_ref = ET.SubElement(listOfProducts, "speciesReference", {
                            "species": self.species[comp_id]["sbml_id"],
                            "stoichiometry": str(stoich),
                            "constant": "true"
                        })
            
            # Add kinetic law (default mass action)
            kineticLaw = ET.SubElement(reaction, "kineticLaw")
            math = ET.SubElement(kineticLaw, "math", {"xmlns": "http://www.w3.org/1998/Math/MathML"})
            
            # Default to mass action kinetics
            if reaction_data.get("reversible", False):
                # Reversible mass action
                ET.SubElement(math, "apply").text = "(k_forward * PRODUCT(reactants) - k_reverse * PRODUCT(products))"
            else:
                # Irreversible mass action
                ET.SubElement(math, "apply").text = "k_forward * PRODUCT(reactants)"
                
            # Add parameters for reaction rates
            listOfParameters = ET.SubElement(kineticLaw, "listOfParameters")
            ET.SubElement(listOfParameters, "parameter", {
                "id": f"k_forward_{sbml_id}",
                "value": "1.0",
                "constant": "true"
            })
            
            if reaction_data.get("reversible", False):
                ET.SubElement(listOfParameters, "parameter", {
                    "id": f"k_reverse_{sbml_id}",
                    "value": "0.1",
                    "constant": "true"
                })
                
        return sbml
    
    def save_sbml(self, output_file: str = None) -> str:
        """
        Convert and save the KEGG pathway to SBML format.
        
        Args:
            output_file: Path to save the SBML file (default: pathway_id.xml)
            
        Returns:
            Path to the saved SBML file
        """
        self.parse_pathway()
        sbml_model = self.create_sbml_model()
        
        # Convert ElementTree to string with pretty formatting
        xml_str = ET.tostring(sbml_model, encoding='utf-8')
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ")
        
        # Save to file
        if output_file is None:
            output_file = f"{self.pathway_id}.xml"
            
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)
            
        return output_file


def main():
    """Main function to run the converter from command line."""
    parser = argparse.ArgumentParser(description='Convert KEGG pathway to SBML format')
    parser.add_argument('-p', '--pathway', default='map00790',
                        help='KEGG Pathway ID (default: map00790 - Folate Biosynthesis)')
    parser.add_argument('-o', '--output', default=None,
                        help='Output file name (default: pathway_id.xml)')
    args = parser.parse_args()

    converter = KEGGtoSBMLConverter(args.pathway)
    output_file = converter.save_sbml(args.output)
    
    print(f"Conversion complete. SBML file saved as: {output_file}")


if __name__ == "__main__":
    main()
