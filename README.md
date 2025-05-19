# KGML to SBML Converter

This Python script converts pathway data from KGML (KEGG Markup Language) format to SBML (Systems Biology Markup Language) format. It is designed to parse a local KGML file, extract pathway information including compounds, genes/enzymes, and reactions, and then represent this information in a structured SBML L3V2 file.

## How it Works

The converter performs the following steps:

1.  **Parses KGML:** Reads the input `.kgml` file using Python's `xml.etree.ElementTree`.
    *   Extracts overall pathway information (ID, name, organism) from the `<pathway>` element.
2.  **Identifies Species:**
    *   Processes `<entry type="compound">` elements to create SBML species for metabolites.
    *   Processes `<entry type="gene">` and `<entry type="ortholog">` elements to create SBML species representing enzymes or catalytic functions. These are typically used as modifiers in reactions.
    *   Unique SBML IDs are generated for all species.
3.  **Defines Reactions:**
    *   Processes `<reaction>` elements in the KGML file.
    *   Identifies substrates and products using their KEGG compound IDs (from the `name` attribute of `<substrate>`/`<product>` sub-elements). Stoichiometry is assumed to be 1 for all participants.
    *   Determines reaction reversibility based on the `type` attribute of the `<reaction>` element.
    *   Associates enzymes (derived from `gene` or `ortholog` entries linked via their `reaction` attribute) as SBML modifiers to the corresponding reactions.
4.  **Generates SBML:**
    *   Constructs an SBML Level 3 Version 2 document.
    *   Creates a single default compartment where all species and reactions reside.
    *   Populates `listOfSpecies` with all identified compounds and enzymes.
    *   Populates `listOfReactions` with details for each reaction, including reactants, products, and enzyme modifiers.
    *   Adds basic MIRIAM-compliant annotations (using `bqbiol:is` or `bqbiol:isDescribedBy`) to link SBML elements back to their KEGG identifiers (e.g., `kegg.compound`, `kegg.reaction`, `kegg.genes`, `kegg.orthology`).
    *   Includes a default mass-action kinetic law for each reaction. Parameters `kf` (forward rate constant) and `kr` (reverse rate constant, if reversible) are added with placeholder values.
5.  **Saves Output:** The resulting SBML model is written to the specified output file with pretty-printed XML formatting.

## Requirements

*   Python 3.6+
*   Standard Python libraries: `xml.etree.ElementTree`, `xml.dom.minidom`, `argparse`, `datetime`, `re`. No external packages are strictly required for the core conversion functionality from a local KGML file.

## How to Use

Run the script from the command line:

```bash
python kgml_to_sbml_converter.py -i path/to/your/input.kgml -o path/to/your/output.xml
Use code with caution.
Markdown
Arguments:

-i, --input FILE_PATH: Required. Path to the input KGML file.
-o, --output FILE_PATH: Required. Path where the output SBML file will be saved (e.g., model.xml).
Example:

python kgml_to_sbml_converter.py -i hsa00010.kgml -o hsa00010_sbml.xml
Use code with caution.
Bash
This will read hsa00010.kgml and produce an SBML file named hsa00010_sbml.xml.

Input KGML Assumptions
The converter expects a KGML file structure similar to those provided by KEGG, including:

A root <pathway> element with attributes like name, org, title.
<entry> elements with id, name, type (e.g., "compound", "gene", "ortholog"), and a nested <graphics> element for display names.
For type="gene" or type="ortholog", an optional reaction attribute listing associated KEGG Reaction IDs (e.g., rn:RXXXXX).
<reaction> elements with id (internal KGML ID), name (KEGG Reaction ID, e.g., rn:RXXXXX), and type ("reversible" or "irreversible").
Nested <substrate> and <product> elements, each with a name attribute specifying the KEGG Compound ID (e.g., cpd:CXXXXX).
Output SBML Features
SBML Level 3 Version 2.
Single Default Compartment: All species and reactions are placed in a compartment named "default_compartment".
Species: Compounds and enzymes are defined as species.
Reactions: Include reactants, products, and enzyme modifiers.
Mass-Action Kinetics: Default kinetic laws (e.g., kf * [S1] * [S2] - kr * [P1] * [P2]) are generated with placeholder rate constants (kf=0.1, kr=0.01).
MIRIAM Annotations: SBML elements are annotated with URIs pointing to identifiers.org for their corresponding KEGG entries, facilitating interoperability.
Notes: The SBML model includes a note about its generation from the KGML file.
