#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import matplotlib.pyplot as plt
from kegg_to_sbml_converter import KEGGtoSBMLConverter

try:
    import basico as bc
except ImportError:
    print("COPASI basico package not found. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "copasi-basico"])
    import basico as bc


def convert_pathway_to_sbml():
    """Convert KEGG folate biosynthesis pathway to SBML"""
    pathway_id = "map00790"  # Folate biosynthesis
    converter = KEGGtoSBMLConverter(pathway_id)
    sbml_file = converter.save_sbml()
    print(f"Converted {pathway_id} to SBML: {sbml_file}")
    return sbml_file


def load_sbml_in_copasi(sbml_file):
    """Load SBML file into COPASI"""
    # Load the model
    model = bc.load_model(sbml_file)
    print(f"Model '{bc.get_model_name()}' loaded successfully")
    
    # Print model information
    print("\nModel Structure:")
    print("-" * 40)
    print(f"Compartments: {len(bc.get_compartments())}")
    print(f"Species: {len(bc.get_species())}")
    print(f"Reactions: {len(bc.get_reactions())}")
    print(f"Parameters: {len(bc.get_parameters())}")
    
    return model


def visualize_model_network():
    """Visualize the reaction network"""
    print("\nGenerating network visualization...")
    # Get connections between species
    rxns = bc.get_reactions()
    connections = []
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Simple visualization - not as sophisticated as COPASI's built-in one
    # but demonstrates programmatic access to model elements
    
    # Get positions for species (simple circular layout)
    species = bc.get_species()
    num_species = len(species)
    species_positions = {}
    
    import math
    for i, (sid, sinfo) in enumerate(species.iterrows()):
        angle = 2 * math.pi * i / num_species
        r = 5  # radius
        x = r * math.cos(angle)
        y = r * math.sin(angle)
        species_positions[sid] = (x, y)
        
        # Draw species
        name = sinfo['name'] if not pd.isna(sinfo['name']) else sid
        name = name[:20] + '...' if len(name) > 20 else name
        ax.plot(x, y, 'o', markersize=10)
        ax.text(x, y+0.3, name, ha='center', va='bottom', fontsize=8)
    
    # Draw reactions
    for rid, rinfo in rxns.iterrows():
        # Get reactants and products
        reactants = bc.get_reaction_parameters(rid, include_species=True)
        reactants = [r for r in reactants.index if '_reference_' in r]
        
        substrates = []
        products = []
        
        for ref in reactants:
            species_id = ref.split('_reference_')[1]
            role = reactants.loc[ref, 'role']
            
            if role == 'substrate' and species_id in species_positions:
                substrates.append(species_id)
            elif role == 'product' and species_id in species_positions:
                products.append(species_id)
        
        # Draw arrows between substrates and products
        for s in substrates:
            for p in products:
                sx, sy = species_positions[s]
                px, py = species_positions[p]
                
                # Draw reaction as midpoint
                mx, my = (sx + px)/2, (sy + py)/2
                
                # Draw lines
                ax.arrow(sx, sy, (mx-sx)*0.8, (my-sy)*0.8, head_width=0.2, 
                         head_length=0.2, fc='blue', ec='blue', alpha=0.5)
                ax.arrow(mx, my, (px-mx)*0.8, (py-my)*0.8, head_width=0.2, 
                         head_length=0.2, fc='green', ec='green', alpha=0.5)
                
                # Add reaction name
                name = rinfo['name'] if not pd.isna(rinfo['name']) else rid
                name = name[:20] + '...' if len(name) > 20 else name
                ax.text(mx, my, name, ha='center', va='center', fontsize=6, 
                        bbox=dict(facecolor='white', alpha=0.7, pad=2))
    
    ax.set_title(f"Reaction Network: {bc.get_model_name()}")
    ax.axis('equal')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig('model_network.png', dpi=300)
    plt.show()
    
    print("Network visualization saved as 'model_network.png'")


def run_time_course():
    """Run a time course simulation"""
    print("\nRunning time course simulation...")
    # Set up simulation parameters
    bc.set_time_course(duration=100, intervals=100, output_event=True)
    
    # Run the simulation
    tc_result = bc.run_time_course()
    
    # Plot the results
    if tc_result is not None and not tc_result.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot only a subset of species if there are many
        species = bc.get_species()
        
        # Select up to 10 species to display
        species_to_plot = species.index[:10] if len(species) > 10 else species.index
        
        for s in species_to_plot:
            if s in tc_result.columns:
                ax.plot(tc_result['Time'], tc_result[s], label=s)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Concentration')
        ax.set_title('Time Course Simulation')
        ax.legend(loc='best', fontsize='small')
        
        plt.tight_layout()
        plt.savefig('time_course.png', dpi=300)
        plt.show()
        
        print("Time course simulation completed and saved as 'time_course.png'")
    else:
        print("Time course simulation failed or returned no data")


def main():
    """Main function to demonstrate the workflow"""
    print("KEGG to SBML converter + COPASI example\n")
    
    # 1. Convert KEGG pathway to SBML
    sbml_file = convert_pathway_to_sbml()
    
    # 2. Load SBML into COPASI
    model = load_sbml_in_copasi(sbml_file)
    
    # 3. Visualize model network
    try:
        import pandas as pd
        visualize_model_network()
    except Exception as e:
        print(f"Could not visualize network: {e}")
    
    # 4. Run time course simulation
    try:
        run_time_course()
    except Exception as e:
        print(f"Could not run time course simulation: {e}")
    
    print("\nWorkflow completed successfully!")


if __name__ == "__main__":
    main()
