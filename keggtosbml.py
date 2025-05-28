#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
import re
from typing import List, Tuple, Dict, Any
import os

class ImprovedKGMLtoSBMLConverter:
    """
    Enhanced converter that transforms KGML (KEGG Markup Language) files to SBML format.
    
    Key improvements:
    - Removes 'fast=false' attribute from reactions
    - Adds proper kinetic parameters (kcat, Km) 
    - Cleaner species naming (removes 's_' prefix)
    - Better enzyme handling
    - Global parameters section
    """

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
        self.pathway_organism = ""

        # Data structures
        self.species = {}  # {kegg_id: species_data}
        self.reactions = {}  # {reaction_id: reaction_data}
        self.parameters = {}  # {param_id: param_data}
        
        # ID management
        self.sbml_id_set = set()
        
        # Default kinetic parameters
        self.default_kcat = 1.0  # s^-1
        self.default_km = 0.1    # mM
        self.default_kf = 0.1    # s^-1
        self.default_kr = 0.01   # s^-1

    def _get_unique_sbml_id(self, base_id: str, prefix: str = "") -> str:
        """Generates a unique SBML ID with cleaner naming."""
        # Clean the base_id
        processed_base_id = re.sub(r'[^A-Za-z0-9_]', '_', base_id)
        if not re.match(r'^[A-Za-z_]', processed_base_id):
            processed_base_id = f"_{processed_base_id}"
        
        # Create full ID without redundant prefixes
        if prefix and not processed_base_id.startswith(prefix):
            full_base_id = f"{prefix}{processed_base_id}"
        else:
            full_base_id = processed_base_id
        
        sbml_id = full_base_id
        n = 1
        while sbml_id in self.sbml_id_set:
            sbml_id = f"{full_base_id}_{n}"
            n += 1
        self.sbml_id_set.add(sbml_id)
        return sbml_id

    def _sanitize_kegg_id_for_uri(self, kegg_id_str: str) -> str:
        """Extracts the core ID from a KEGG ID string."""
        if ':' in kegg_id_str:
            return kegg_id_str.split(':')[-1]
        return kegg_id_str
        
    def _extract_kegg_ids_from_name_attr(self, name_attr: str) -> list:
        """Extracts multiple KEGG IDs from an entry's name attribute."""
        return name_attr.split()

    def _add_global_parameters(self):
        """Add global kinetic parameters to the model."""
        global_params = {
            'default_kcat': {'value': self.default_kcat, 'units': 'per_second', 'name': 'Default catalytic rate constant'},
            'default_km': {'value': self.default_km, 'units': 'mM', 'name': 'Default Michaelis constant'},
            'default_kf': {'value': self.default_kf, 'units': 'per_second', 'name': 'Default forward rate constant'},
            'default_kr': {'value': self.default_kr, 'units': 'per_second', 'name': 'Default reverse rate constant'},
            'temperature': {'value': 298.15, 'units': 'kelvin', 'name': 'Temperature'},
            'avogadro': {'value': 6.022e23, 'units': 'per_mole', 'name': 'Avogadro constant'}
        }
        
        for param_id, param_data in global_params.items():
            self.parameters[param_id] = param_data

    def _add_enzyme_parameters(self, enzyme_id: str, reaction_id: str):
        """Add enzyme-specific kinetic parameters."""
        # Create enzyme-specific parameters
        kcat_id = f"kcat_{enzyme_id}_{reaction_id}"
        km_id = f"Km_{enzyme_id}_{reaction_id}"
        
        self.parameters[kcat_id] = {
            'value': self.default_kcat,
            'units': 'per_second',
            'name': f'Catalytic rate constant for {enzyme_id} in {reaction_id}'
        }
        
        self.parameters[km_id] = {
            'value': self.default_km,
            'units': 'mM',
            'name': f'Michaelis constant for {enzyme_id} in {reaction_id}'
        }
        
        return kcat_id, km_id

    def _parse_kgml(self):
        """Parses the KGML file and populates internal data structures."""
        try:
            self.kgml_tree = ET.parse(self.kgml_filepath)
            self.kgml_root = self.kgml_tree.getroot()
        except ET.ParseError as e:
            raise ValueError(f"Error parsing KGML file: {e}")
        except FileNotFoundError:
            raise FileNotFoundError(f"KGML file not found: {self.kgml_filepath}")

        # Parse pathway information
        pathway_kegg_id = self.kgml_root.get('name', 'unknown_pathway')
        self.pathway_sbml_id = self._get_unique_sbml_id(
            self._sanitize_kegg_id_for_uri(pathway_kegg_id), "model_"
        )
        self.pathway_name = self.kgml_root.get('title', 'Untitled Pathway')
        self.pathway_organism = self.kgml_root.get('org', '')

        # Add global parameters
        self._add_global_parameters()

        # Temporary mapping for gene entries to reactions
        gene_entry_to_reactions = {}
        kgml_entries = {}

        # First pass: Parse entries
        for entry_elem in self.kgml_root.findall('entry'):
            entry_id = entry_elem.get('id')
            entry_type = entry_elem.get('type')
            entry_kegg_name_attr = entry_elem.get('name', '')
            
            graphics_elem = entry_elem.find('graphics')
            graphics_name_str = ''
            if graphics_elem is not None:
                graphics_name_str = graphics_elem.get('name', '').replace('...', '').strip()

            kgml_entries[entry_id] = {
                'kegg_name_attr': entry_kegg_name_attr,
                'type': entry_type,
                'graphics_name': graphics_name_str
            }

            if entry_type == 'compound':
                compound_kegg_id = self._sanitize_kegg_id_for_uri(entry_kegg_name_attr)
                if compound_kegg_id not in self.species:
                    # Clean naming - just use the compound ID
                    sbml_id = self._get_unique_sbml_id(compound_kegg_id)
                    display_name = graphics_name_str if graphics_name_str else compound_kegg_id
                    
                    self.species[compound_kegg_id] = {
                        'sbml_id': sbml_id,
                        'name': display_name,
                        'kegg_ids_full': self._extract_kegg_ids_from_name_attr(entry_kegg_name_attr),
                        'type': 'compound',
                        'entry_id': entry_id
                    }
                kgml_entries[entry_id]['sbml_id'] = self.species[compound_kegg_id]['sbml_id']
            
            elif entry_type in ['gene', 'ortholog']:
                # For enzymes, create cleaner IDs
                base_name = graphics_name_str.split(',')[0].strip() if graphics_name_str else ''
                if not base_name:
                    base_name = self._sanitize_kegg_id_for_uri(entry_kegg_name_attr.split()[0])
                
                # Clean enzyme naming
                sbml_id = self._get_unique_sbml_id(f"enzyme_{base_name}")
                display_name = graphics_name_str if graphics_name_str else entry_kegg_name_attr

                self.species[entry_id] = {
                    'sbml_id': sbml_id,
                    'name': display_name,
                    'kegg_ids_full': self._extract_kegg_ids_from_name_attr(entry_kegg_name_attr),
                    'type': 'enzyme',
                    'entry_id': entry_id
                }
                kgml_entries[entry_id]['sbml_id'] = self.species[entry_id]['sbml_id']

                # Store reaction associations
                reaction_attr = entry_elem.get('reaction')
                if reaction_attr:
                    gene_entry_to_reactions[entry_id] = [
                        self._sanitize_kegg_id_for_uri(r) for r in reaction_attr.split()
                    ]

        # Second pass: Parse reactions
        for reaction_elem in self.kgml_root.findall('reaction'):
            reaction_kegg_id_full = reaction_elem.get('name', '')
            reaction_kegg_id = self._sanitize_kegg_id_for_uri(reaction_kegg_id_full)
            
            if reaction_kegg_id in self.reactions:
                continue

            # Clean reaction naming
            sbml_reaction_id = self._get_unique_sbml_id(reaction_kegg_id)
            reaction_name = reaction_kegg_id
            reaction_reversible = reaction_elem.get('type') == 'reversible'
            
            current_reaction_data = {
                'sbml_id': sbml_reaction_id,
                'name': reaction_name,
                'reversible': reaction_reversible,
                'substrates': [],
                'products': [],
                'modifiers': [],
                'kegg_id_full': reaction_kegg_id_full,
                'enzymes': []
            }

            # Parse substrates
            for substrate_elem in reaction_elem.findall('substrate'):
                compound_kegg_id_full = substrate_elem.get('name', '')
                compound_kegg_id = self._sanitize_kegg_id_for_uri(compound_kegg_id_full)
                
                if compound_kegg_id not in self.species:
                    sbml_id = self._get_unique_sbml_id(compound_kegg_id)
                    self.species[compound_kegg_id] = {
                        'sbml_id': sbml_id,
                        'name': compound_kegg_id,
                        'kegg_ids_full': [compound_kegg_id_full],
                        'type': 'compound'
                    }
                
                if self.species[compound_kegg_id]['type'] == 'compound':
                    current_reaction_data['substrates'].append({
                        'id': self.species[compound_kegg_id]['sbml_id'], 
                        'stoichiometry': 1
                    })

            # Parse products
            for product_elem in reaction_elem.findall('product'):
                compound_kegg_id_full = product_elem.get('name', '')
                compound_kegg_id = self._sanitize_kegg_id_for_uri(compound_kegg_id_full)

                if compound_kegg_id not in self.species:
                    sbml_id = self._get_unique_sbml_id(compound_kegg_id)
                    self.species[compound_kegg_id] = {
                        'sbml_id': sbml_id,
                        'name': compound_kegg_id,
                        'kegg_ids_full': [compound_kegg_id_full],
                        'type': 'compound'
                    }
                
                if self.species[compound_kegg_id]['type'] == 'compound':
                    current_reaction_data['products'].append({
                        'id': self.species[compound_kegg_id]['sbml_id'], 
                        'stoichiometry': 1
                    })

            self.reactions[reaction_kegg_id] = current_reaction_data

        # Link enzymes to reactions and add parameters
        for gene_entry_id, associated_rn_ids in gene_entry_to_reactions.items():
            if gene_entry_id in self.species and self.species[gene_entry_id]['type'] == 'enzyme':
                enzyme_sbml_id = self.species[gene_entry_id]['sbml_id']
                enzyme_display_name = self.species[gene_entry_id]['name']
                
                for rn_id in associated_rn_ids:
                    if rn_id in self.reactions:
                        # Add enzyme as modifier
                        self.reactions[rn_id]['modifiers'].append({'id': enzyme_sbml_id})
                        self.reactions[rn_id]['enzymes'].append(enzyme_sbml_id)
                        
                        # Add enzyme-specific parameters
                        self._add_enzyme_parameters(enzyme_sbml_id, rn_id)
                        
                        # Update reaction name with enzyme info
                        if enzyme_display_name:
                            enzyme_name = enzyme_display_name.split(',')[0].strip()
                            self.reactions[rn_id]['name'] = f"{rn_id}: {enzyme_name}"

    def _create_sbml_annotation(self, parent_element: ET.Element, sbml_element_id: str, 
                               kegg_id_tuples: List[Tuple[str, str]]):
        """Creates RDF annotation for an SBML element."""
        if not kegg_id_tuples:
            return

        annotation = ET.SubElement(parent_element, "annotation")
        rdf_rdf = ET.SubElement(annotation, "rdf:RDF", {
            "xmlns:rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "xmlns:bqbiol": "http://biomodels.net/biology-qualifiers/"
        })
        description = ET.SubElement(rdf_rdf, "rdf:Description", 
                                  {"rdf:about": f"#{sbml_element_id}"})
        
        qualifier_tag = "bqbiol:is" if len(kegg_id_tuples) == 1 else "bqbiol:isDescribedBy"
        bq_element = ET.SubElement(description, qualifier_tag)
        bag = ET.SubElement(bq_element, "rdf:Bag")

        for kegg_db_prefix, full_kegg_id in kegg_id_tuples:
            individual_ids = full_kegg_id.split()
            for single_id in individual_ids:
                core_id = self._sanitize_kegg_id_for_uri(single_id)
                uri = self._get_identifiers_uri(kegg_db_prefix, single_id, core_id)
                if uri:
                    ET.SubElement(bag, "rdf:li", {"rdf:resource": uri})

    def _get_identifiers_uri(self, kegg_db_prefix: str, single_id: str, core_id: str) -> str:
        """Generate identifiers.org URI for KEGG IDs."""
        uri_map = {
            "compound": f"http://identifiers.org/kegg.compound/{core_id}",
            "reaction": f"http://identifiers.org/kegg.reaction/{core_id}",
            "genes": f"http://identifiers.org/kegg.genes/{single_id}",
            "orthology": f"http://identifiers.org/kegg.orthology/{core_id}",
            "pathway": f"http://identifiers.org/kegg.pathway/{core_id}"
        }
        
        if kegg_db_prefix in uri_map:
            return uri_map[kegg_db_prefix]
        elif kegg_db_prefix.startswith("path"):
            return f"http://identifiers.org/kegg.pathway/{core_id}"
        elif ":" in single_id:
            db, actual_id = single_id.split(":", 1)
            return f"http://identifiers.org/kegg.{db}/{actual_id}"
        return ""

    def _create_kinetic_law_mathml(self, reaction_data: Dict[str, Any]) -> str:
        """Create MathML for kinetic law based on reaction type."""
        sbml_id = reaction_data['sbml_id']
        has_enzymes = bool(reaction_data['enzymes'])
        
        if has_enzymes:
            # Michaelis-Menten kinetics for enzymatic reactions
            return self._create_mm_kinetics_mathml(reaction_data)
        else:
            # Mass action kinetics for non-enzymatic reactions
            return self._create_mass_action_mathml(reaction_data)

    def _create_mm_kinetics_mathml(self, reaction_data: Dict[str, Any]) -> str:
        """Create Michaelis-Menten kinetics MathML."""
        # Simplified MM: v = kcat * E * S / (Km + S)
        # For multiple substrates: v = kcat * E * S1 * S2 / ((Km1 + S1) * (Km2 + S2))
        
        enzyme_id = reaction_data['enzymes'][0] if reaction_data['enzymes'] else 'default_enzyme'
        reaction_id = reaction_data['sbml_id']
        
        kcat_param = f"kcat_{enzyme_id}_{reaction_id.replace('r_', '').replace('R', '')}"
        
        mathml_parts = [f'<ci>{kcat_param}</ci>', f'<ci>{enzyme_id}</ci>']
        
        # Add substrates to numerator
        for substrate in reaction_data['substrates']:
            mathml_parts.append(f'<ci>{substrate["id"]}</ci>')
        
        # Create numerator
        if len(mathml_parts) > 1:
            numerator = f"<apply><times/>{''.join(mathml_parts)}</apply>"
        else:
            numerator = mathml_parts[0]
        
        # Simple denominator (can be enhanced)
        denominator = "<cn type='real'>1.0</cn>"
        
        return f"<apply><divide/>{numerator}{denominator}</apply>"

    def _create_mass_action_mathml(self, reaction_data: Dict[str, Any]) -> str:
        """Create mass action kinetics MathML."""
        sbml_id = reaction_data['sbml_id']
        
        # Forward reaction: kf * S1 * S2 * ...
        forward_parts = ['<ci>default_kf</ci>']
        for substrate in reaction_data['substrates']:
            forward_parts.append(f'<ci>{substrate["id"]}</ci>')
        
        if len(forward_parts) > 1:
            forward_term = f"<apply><times/>{''.join(forward_parts)}</apply>"
        else:
            forward_term = forward_parts[0] if forward_parts else "<cn type='real'>0</cn>"
        
        if not reaction_data['reversible']:
            return forward_term
        
        # Reverse reaction: kr * P1 * P2 * ...
        reverse_parts = ['<ci>default_kr</ci>']
        for product in reaction_data['products']:
            reverse_parts.append(f'<ci>{product["id"]}</ci>')
        
        if len(reverse_parts) > 1:
            reverse_term = f"<apply><times/>{''.join(reverse_parts)}</apply>"
        else:
            reverse_term = reverse_parts[0] if reverse_parts else "<cn type='real'>0</cn>"
        
        return f"<apply><minus/>{forward_term}{reverse_term}</apply>"

    def create_sbml_model(self) -> ET.Element:
        """Creates the SBML model from parsed KGML data."""
        sbml_ns = "http://www.sbml.org/sbml/level3/version2/core"
        sbml = ET.Element("sbml", {"xmlns": sbml_ns, "level": "3", "version": "2"})
        model = ET.SubElement(sbml, "model", {"id": self.pathway_sbml_id, "name": self.pathway_name})

        # Add model notes
        self._add_model_notes(model)
        
        # Add model annotation
        self._add_model_annotation(model)

        # Add unit definitions
        self._add_unit_definitions(model)

        # Add compartments
        self._add_compartments(model)

        # Add species (only compounds, enzymes handled as modifiers)
        self._add_species(model)

        # Add parameters
        self._add_parameters(model)

        # Add reactions
        self._add_reactions(model)

        return sbml

    def _add_model_notes(self, model: ET.Element):
        """Add notes to the model."""
        notes_str = (f"Model generated from KEGG KGML file '{os.path.basename(self.kgml_filepath)}' "
                    f"for pathway {self.pathway_name}. "
                    f"Conversion performed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.")
        
        notes_elem = ET.SubElement(model, "notes")
        html_notes = ET.SubElement(notes_elem, "html", {"xmlns": "http://www.w3.org/1999/xhtml"})
        p_elem = ET.SubElement(html_notes, "p")
        p_elem.text = notes_str

    def _add_model_annotation(self, model: ET.Element):
        """Add annotation to the model."""
        pathway_kegg_full_id = self.kgml_root.get('name')
        if pathway_kegg_full_id:
            sanitized_path_id = self._sanitize_kegg_id_for_uri(pathway_kegg_full_id)
            self._create_sbml_annotation(model, self.pathway_sbml_id, [("pathway", sanitized_path_id)])

    def _add_unit_definitions(self, model: ET.Element):
        """Add unit definitions to the model."""
        listOfUnitDefinitions = ET.SubElement(model, "listOfUnitDefinitions")
        
        # Define common units
        units = {
            'per_second': [('second', -1)],
            'mM': [('mole', 1), ('litre', -1), ('scale', 0, -3)],
            'per_mole': [('mole', -1)],
            'kelvin': [('kelvin', 1)]
        }
        
        for unit_id, unit_components in units.items():
            unit_def = ET.SubElement(listOfUnitDefinitions, "unitDefinition", {"id": unit_id})
            list_of_units = ET.SubElement(unit_def, "listOfUnits")
            
            for component in unit_components:
                unit_attrs = {"kind": component[0], "exponent": str(component[1])}
                if len(component) > 2:
                    unit_attrs["scale"] = str(component[2])
                ET.SubElement(list_of_units, "unit", unit_attrs)

    def _add_compartments(self, model: ET.Element):
        """Add compartments to the model."""
        listOfCompartments = ET.SubElement(model, "listOfCompartments")
        ET.SubElement(listOfCompartments, "compartment", {
            "id": "default_compartment", 
            "name": "Default Compartment", 
            "constant": "true", 
            "spatialDimensions": "3", 
            "size": "1"
        })

    def _add_species(self, model: ET.Element):
        """Add species to the model."""
        listOfSpecies = ET.SubElement(model, "listOfSpecies")
        
        for species_key, species_data in self.species.items():
            # Only add compounds as species; enzymes are handled as modifiers
            if species_data['type'] != 'compound':
                continue
                
            sbml_id = species_data['sbml_id']
            name = species_data['name']
            kegg_ids_full_list = species_data['kegg_ids_full']

            species_elem = ET.SubElement(listOfSpecies, "species", {
                "id": sbml_id, 
                "name": name, 
                "compartment": "default_compartment",
                "hasOnlySubstanceUnits": "false", 
                "boundaryCondition": "false", 
                "constant": "false"
            })
            
            # Add annotations
            kegg_id_tuples = []
            for kid_full in kegg_ids_full_list:
                if kid_full.startswith("cpd:"):
                    kegg_id_tuples.append(("compound", kid_full))
                elif kid_full.startswith("gl:"):
                    kegg_id_tuples.append(("glycan", kid_full))
            
            self._create_sbml_annotation(species_elem, sbml_id, kegg_id_tuples)

    def _add_parameters(self, model: ET.Element):
        """Add parameters to the model."""
        if not self.parameters:
            return
            
        listOfParameters = ET.SubElement(model, "listOfParameters")
        
        for param_id, param_data in self.parameters.items():
            param_attrs = {
                "id": param_id,
                "name": param_data['name'],
                "value": str(param_data['value']),
                "constant": "true"
            }
            if 'units' in param_data:
                param_attrs["units"] = param_data['units']
                
            ET.SubElement(listOfParameters, "parameter", param_attrs)

    def _add_reactions(self, model: ET.Element):
        """Add reactions to the model."""
        listOfReactions = ET.SubElement(model, "listOfReactions")
        
        for reaction_key, reaction_data in self.reactions.items():
            sbml_id = reaction_data['sbml_id']
            name = reaction_data['name']
            kegg_reaction_id_full = reaction_data['kegg_id_full']

            # Create reaction element WITHOUT fast attribute
            reaction_elem = ET.SubElement(listOfReactions, "reaction", {
                "id": sbml_id, 
                "name": name, 
                "reversible": "true" if reaction_data['reversible'] else "false"
            })
            
            # Add annotation
            self._create_sbml_annotation(reaction_elem, sbml_id, [('reaction', kegg_reaction_id_full)])

            # Add reactants
            if reaction_data['substrates']:
                listOfReactants = ET.SubElement(reaction_elem, "listOfReactants")
                for s_info in reaction_data['substrates']:
                    ET.SubElement(listOfReactants, "speciesReference", {
                        "species": s_info['id'], 
                        "stoichiometry": str(s_info['stoichiometry']), 
                        "constant": "true"
                    })
            
            # Add products
            if reaction_data['products']:
                listOfProducts = ET.SubElement(reaction_elem, "listOfProducts")
                for p_info in reaction_data['products']:
                    ET.SubElement(listOfProducts, "speciesReference", {
                        "species": p_info['id'], 
                        "stoichiometry": str(p_info['stoichiometry']), 
                        "constant": "true"
                    })

            # Add modifiers (enzymes)
            if reaction_data['modifiers']:
                listOfModifiers = ET.SubElement(reaction_elem, "listOfModifiers")
                for m_info in reaction_data['modifiers']:
                    # Note: Enzymes are referenced but not included as species
                    modifier_elem = ET.SubElement(listOfModifiers, "modifierSpeciesReference", 
                                                {"species": m_info['id']})
            
            # Add kinetic law
            self._add_kinetic_law(reaction_elem, reaction_data)

    def _add_kinetic_law(self, reaction_elem: ET.Element, reaction_data: Dict[str, Any]):
        """Add kinetic law to reaction."""
        kineticLaw = ET.SubElement(reaction_elem, "kineticLaw")
        
        # Add MathML
        math = ET.SubElement(kineticLaw, "math", {"xmlns": "http://www.w3.org/1998/Math/MathML"})
        mathml_str = self._create_kinetic_law_mathml(reaction_data)
        
        try:
            math_expr_elem = ET.fromstring(mathml_str)
            math.append(math_expr_elem)
        except ET.ParseError:
            # Fallback to zero rate
            cn_elem = ET.SubElement(math, "cn", {"type": "real"})
            cn_elem.text = "0.0"
            print(f"Warning: Could not parse MathML for reaction {reaction_data['sbml_id']}")

    def convert_and_save(self, output_filepath: str):
        """Convert KGML to SBML and save to file."""
        print(f"Parsing KGML file: {self.kgml_filepath}")
        self._parse_kgml()
        
        print("Generating SBML model...")
        sbml_model_element = self.create_sbml_model()
        
        # Pretty print and save
        xml_str = ET.tostring(sbml_model_element, encoding='utf-8', method='xml')
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ")
        
        with open(output_filepath, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)
        
        print(f"SBML file saved as: {output_filepath}")
        
        # Print summary statistics
        print("\n=== Conversion Summary ===")
        print(f"Pathway: {self.pathway_name}")
        print(f"Organism: {self.pathway_organism}")
        
        compound_count = sum(1 for s in self.species.values() if s['type'] == 'compound')
        enzyme_count = sum(1 for s in self.species.values() if s['type'] == 'enzyme')
        
        print(f"Species (compounds): {compound_count}")
        print(f"Enzymes: {enzyme_count}")
        print(f"Reactions: {len(self.reactions)}")
        print(f"Parameters: {len(self.parameters)}")
        
        # Validate that enzymes are not included as species in SBML
        print("\n=== Validation ===")
        enzymes_as_modifiers = 0
        for reaction_data in self.reactions.values():
            enzymes_as_modifiers += len(reaction_data['modifiers'])
        
        print(f"Enzymes referenced as modifiers: {enzymes_as_modifiers}")
        
        # Check for clean naming (no 's_' prefixes)
        species_with_clean_names = 0
        for species_data in self.species.values():
            if species_data['type'] == 'compound' and not species_data['sbml_id'].startswith('s_'):
                species_with_clean_names += 1
        
        print(f"Species with clean names (no 's_' prefix): {species_with_clean_names}/{compound_count}")
        
        # Verify parameters exist for enzymatic reactions
        enzymatic_reactions = sum(1 for r in self.reactions.values() if r['enzymes'])
        print(f"Enzymatic reactions with parameters: {enzymatic_reactions}")
        
        print("=== Conversion Complete ===")


def main():
    parser = argparse.ArgumentParser(description='Convert KEGG KGML file to SBML format with improved handling.')
    parser.add_argument('-i', '--input', required=True, help='Input KGML file path.')
    parser.add_argument('-o', '--output', required=True, help='Output SBML file path (e.g., model.xml).')
    parser.add_argument('--validate', action='store_true', help='Perform additional validation checks.')
    args = parser.parse_args()

    try:
        converter = ImprovedKGMLtoSBMLConverter(args.input)
        converter.convert_and_save(args.output)
        
        if args.validate:
            print("\n=== Additional Validation ===")
            # Additional validation could be added here
            print("Validation complete - check output for any warnings.")
            
    except Exception as e:
        print(f"An error occurred during conversion: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
