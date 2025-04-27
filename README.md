# KEGG to SBML Converter

A comprehensive Python tool for parsing KEGG pathways and converting them to SBML (Systems Biology Markup Language) format for use in systems biology simulation platforms like VCell and COPASI.

## Overview

KEGG to SBML Converter bridges the gap between KEGG (Kyoto Encyclopedia of Genes and Genomes) pathway databases and computational systems biology modeling tools. It automatically extracts metabolic pathway information from KEGG and reformats it into the widely-used SBML standard, enabling simulation, visualization, and analysis in platforms like VCell and COPASI.

## Features

- **Direct KEGG API Access**: Fetches pathway data directly from KEGG's REST API
- **Comprehensive Parsing**: Extracts compounds, reactions, stoichiometry, and reversibility information
- **SBML Level 3 Export**: Creates valid SBML Level 3 Version 1 documents
- **Flexible Pathway Selection**: Works with any KEGG pathway by ID
- **Default Focus on Folate Biosynthesis**: Pre-configured for the Folate Biosynthesis pathway (map00790)
- **Mass Action Kinetics**: Includes default mass action kinetic laws
- **COPASI Integration**: Example scripts for loading models into COPASI
- **Visualization Support**: Helper functions for visualizing the model network
- **Simulation Ready**: Models can be directly simulated in COPASI or VCell

## Requirements

- Python 3.7+
- Required Python packages:
  - requests
  - xml.etree.ElementTree (part of standard library)
  - minidom (part of standard library)
  
- For the example COPASI integration:
  - copasi-basico
  - pandas
  - numpy
  - matplotlib

## Installation

### Direct Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/kegg-to-sbml.git
cd kegg-to-sbml
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

### As a Python Package

```bash
pip install kegg-to-sbml
```

## Usage

### Basic Command Line Usage

Convert a KEGG pathway to SBML format with the default command:

```bash
python kegg_to_sbml_converter.py
```

This will convert the default Folate Biosynthesis pathway (map00790) and save it as "map00790.xml".

#### Specifying a Different Pathway

```bash
python kegg_to_sbml_converter.py --pathway map00010
```

This will convert the Glycolysis pathway (map00010).

#### Specifying the Output File

```bash
python kegg_to_sbml_converter.py --pathway map00790 --output folate_biosynthesis.xml
```

#### Full Command-line Options

```bash
python kegg_to_sbml_converter.py --help
```

Output:
```
usage: kegg_to_sbml_converter.py [-h] [-p PATHWAY] [-o OUTPUT]

Convert KEGG pathway to SBML format

optional arguments:
  -h, --help            show this help message and exit
  -p PATHWAY, --pathway PATHWAY
                        KEGG Pathway ID (default: map00790 - Folate Biosynthesis)
  -o OUTPUT, --output OUTPUT
                        Output file name (default: pathway_id.xml)
```

### Python API Usage

You can also use the converter as a Python module in your own scripts:

```python
from kegg_to_sbml_converter import KEGGtoSBMLConverter

# Create converter for a specific pathway
converter = KEGGtoSBMLConverter("map00790")  # Folate biosynthesis

# Parse and convert to SBML
converter.parse_pathway()
sbml_file = converter.save_sbml("output.xml")

print(f"SBML file saved to: {sbml_file}")
```

### Examples

#### Converting and Analyzing a Pathway with COPASI

The repository includes an example script demonstrating how to:
1. Convert a KEGG pathway to SBML
2. Load it into COPASI using the basico package
3. Visualize the model
4. Run a time course simulation

Run the example script:

```bash
python kegg_to_sbml_example.py
```

This will:
- Convert the Folate Biosynthesis pathway to SBML
- Load the SBML model into COPASI
- Generate a network visualization
- Run a time course simulation
- Save visualization plots as PNG files

## Technical Details

### KEGG Database Access

The tool accesses KEGG data through its REST API:
- Base URL: http://rest.kegg.jp
- Pathway data: /get/{pathway_id}
- Compound data: /get/{compound_id}
- Reaction data: /get/{reaction_id}

No API key is required for basic KEGG API access, but usage should respect KEGG's [terms of use](https://www.kegg.jp/kegg/legal.html).

### Parsing Process

The converter follows these steps to parse KEGG data:

1. **Pathway Retrieval**: Fetches the complete pathway entry
2. **Extract Compounds**: Parses the COMPOUND section to identify all metabolites
3. **Extract Reactions**: Parses the REACTION section to identify all reactions
4. **Compound Details**: For each compound, retrieves detailed information including name, formula, and charge
5. **Reaction Details**: For each reaction, retrieves:
   - Name and description
   - Complete reaction equation
   - Substrates and products with stoichiometry
   - Reversibility information

### SBML Generation

The SBML file generation includes:

1. **SBML Level 3 Structure**: Creates a valid SBML Level 3 Version 1 document
2. **Model Element**: Creates the main model with the pathway name
3. **Compartments**: Creates a default "cytosol" compartment
4. **Species**: Adds all compounds as species with:
   - Unique IDs
   - Names
   - Formula annotations when available
   - Compartment assignments
5. **Reactions**: Adds all reactions with:
   - Reactants and products with stoichiometry
   - Reversibility flags
   - Default mass action kinetics
   - Default parameter values

## Integration with Simulation Tools

### COPASI Integration

The included example shows how to use the COPASI basico package to:

1. **Load SBML**: Import the generated SBML file
2. **Model Analysis**: Access model components programmatically
3. **Visualization**: Create network diagrams
4. **Simulation**: Run time course simulations
5. **Plot Results**: Visualize simulation outcomes

Example code:

```python
import basico as bc

# Load the SBML model
model = bc.load_model("folate_biosynthesis.xml")

# Print model information
print(f"Model: {bc.get_model_name()}")
print(f"Species: {len(bc.get_species())}")
print(f"Reactions: {len(bc.get_reactions())}")

# Run a simulation
bc.set_time_course(duration=100, intervals=100)
tc_result = bc.run_time_course()

# Plot results
import matplotlib.pyplot as plt
species = bc.get_species().index[:5]  # First 5 species
for s in species:
    plt.plot(tc_result['Time'], tc_result[s], label=s)
plt.legend()
plt.show()
```

### VCell Integration

To import the generated SBML file into VCell:

1. Open VCell (either desktop application or web-based version)
2. Create a new BioModel or MathModel
3. Import the SBML file:
   - Go to File → Import → SBML
   - Select your generated SBML file
4. VCell will import the model components
5. You can then:
   - Add spatial information if needed
   - Configure simulation parameters
   - Run simulations
   - Analyze results

## Customization

### Adding Custom Kinetics

The default implementation uses mass action kinetics. To customize:

1. Modify the `create_sbml_model` method in the `KEGGtoSBMLConverter` class
2. Change the math element generation in the kinetic law section
3. Add custom parameters as needed

Example for Michaelis-Menten kinetics:

```python
# In the kineticLaw section for irreversible reactions:
ET.SubElement(math, "apply").text = "Vmax * S / (Km + S)"

# Add parameters
ET.SubElement(listOfParameters, "parameter", {
    "id": f"Vmax_{sbml_id}",
    "value": "1.0",
    "constant": "true"
})
ET.SubElement(listOfParameters, "parameter", {
    "id": f"Km_{sbml_id}",
    "value": "0.5",
    "constant": "true"
})
```

### Adding Annotations

You can extend the tool to include more biological annotations:

```python
# Add EC number annotation to a reaction
if reaction_data.get("ec_number"):
    annotation = ET.SubElement(reaction, "annotation")
    rdf_desc = ET.SubElement(annotation, "rdf:Description")
    ET.SubElement(rdf_desc, "ec").text = reaction_data["ec_number"]
```

## Limitations

- **Default Kinetics**: Uses simple mass action kinetics as a placeholder
- **Limited Compartmentalization**: Assigns all species to a default compartment
- **No Enzyme Information**: Currently doesn't include detailed enzyme data
- **Simplified MathML**: Uses text representation instead of full MathML structure
- **No Regulation**: Does not capture regulatory effects or modifiers
- **Initial Concentrations**: Uses default initial concentrations for all species



