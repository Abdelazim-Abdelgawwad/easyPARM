#!/bin/bash
###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              # 
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
#   /$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$| $$  \ $$| $$  \ $$| $$  \ $$| $$$$  /$$$$                             #
#  /$$__  $$ |____  $$ /$$_____/| $$  | $$| $$$$$$$/| $$$$$$$$| $$$$$$$/| $$ $$/$$ $$                             #
# | $$$$$$$$  /$$$$$$$|  $$$$$$ | $$  | $$| $$____/ | $$__  $$| $$__  $$| $$  $$$| $$                             #
# | $$_____/ /$$__  $$ \____  $$| $$  | $$| $$      | $$  | $$| $$  \ $$| $$\  $ | $$                             #
# |  $$$$$$$|  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$      | $$  | $$| $$  | $$| $$ \/  | $$                             #
#  \_______/ \_______/|_______/  \____  $$|__/      |__/  |__/|__/  |__/|__/     |__/                             #
#                               /$$  | $$                                                                         #
#                              |  $$$$$$/              Ver. 3.25 - 14 April 2025                                  #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


import re
from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import os

#Update atom types and charges in a MOL2 file based on a charge file.
def update_mol2_file(input_mol2, charge_file, output_mol2=None):
    try:
        # If no output file specified, generate one
        if output_mol2 is None:
            base, ext = os.path.splitext(input_mol2)
            output_mol2 = f"{os.path.basename(input_mol2)}"

        # Read new atom charges and types
        new_data = {}
        with open(charge_file, 'r') as f:
            for index, line in enumerate(f, 1):
                parts = line.split()
                if len(parts) >= 2:
                    charge = float(parts[0])
                    atom_type = parts[1]
                    new_data[index] = {'charge': charge, 'atom_type': atom_type}
        
        # Process MOL2 file
        updated_mol2_lines = []
        with open(input_mol2, 'r') as input_file:
            is_atom_section = False
            for line in input_file:
                # Detect and modify the atom section
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    continue
                elif line.startswith("@<TRIPOS>"):
                    is_atom_section = False
                    updated_mol2_lines.append(line)
                    continue
                
                if is_atom_section and line.strip():
                    # Match MOL2 atom line format using regex
                    atom_match = re.match(
                        r"(\s*\d+\s+)(\S+\s+)(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+\s+)(\d+\s+\S+\s+)(-?\d+\.\d+)",
                        line
                    )
                    if atom_match:
                        atom_id, atom_name, coords, old_atom_type, post_type, old_charge = atom_match.groups()
                        atom_index = int(atom_id.strip())
                        if atom_index in new_data:
                            # Replace atom type and charge
                            new_atom_type = new_data[atom_index]['atom_type']
                            new_charge = new_data[atom_index]['charge']
                            updated_line = f"{atom_id}{atom_name}{coords}{new_atom_type:<11}{post_type}{new_charge:8.6f}\n"
                            updated_mol2_lines.append(updated_line)
                        else:
                            updated_mol2_lines.append(line)
                    else:
                        updated_mol2_lines.append(line)
                else:
                    updated_mol2_lines.append(line)
        
        # Write the updated MOL2 content to the output file
        with open(output_mol2, 'w') as output_file:
            output_file.writelines(updated_mol2_lines)
        
        return output_mol2
    
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        raise
    except PermissionError:
        print("Error: Permission denied when accessing files.")
        raise
    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        raise

#Batch update MOL2 files based on a mapping file.    
def batch_update_mol2_files(charges_mapping_file='charges_all.dat'):
    updated_files = []
    try:
        with open(charges_mapping_file, 'r') as mapping_file:
            for line in mapping_file:
                parts = line.strip().split()
                if len(parts) == 2:
                    input_mol2, charge_file = parts
                    try:
                        updated_file = update_mol2_file(input_mol2, charge_file)
                        updated_files.append(updated_file)
                    except Exception as e:
                        print(f"Failed to update {input_mol2}: {e}")
    
    except FileNotFoundError:
        print(f"Error: Charges mapping file '{charges_mapping_file}' not found.")
    
    return updated_files

#Extracts coordination information for standard protein residues linked to metals
#and their specific peptide bond connections.
def generate_standard_residue_coordination(input_pdb, output_file="coordinated_residues.txt", 
                                           standard_residues=None, 
                                           metals=None, 
                                           distance_cutoff=2.5,
                                           bond_cutoff=1.9):
    # Default standard residues
    if standard_residues is None:
        standard_residues = {
            "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
            "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }
    
    # Default metals to check
    if metals is None:
        metals = [
            'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Na', 'K', 'Li', 'Rb', 'Cs', 'Mg',
            'Ca', 'Sr', 'Ba', 'V', 'Cr', 'Cd', 'Hg', 'Al', 'Ga', 'In', 'Sn', 'Pb',
            'Bi', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
            'Er', 'Tm', 'Yb', 'Lu', 'Fe2', 'Fe3', 'Fe4', 'Cu1', 'Cu2', 'Mn2', 'Mn3',
            'Mn4', 'Co2', 'Co3', 'Ni2', 'Ni3', 'V2', 'V3', 'V4', 'V5'
        ] 

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    coordination_data = []
    peptide_bond_data = []
    nonstandard_standard_bond_data = []
    
    # Sets to store residue numbers for tracking
    metal_coordinated_residues = set()
    standard_residues_to_process = set()
    
    # Check if PDB was loaded correctly
    model_count = len(list(structure.get_models()))
    
    metal_atoms_found = []
    nonstandard_residues_found = set()
    standard_residues_found = set()
    
    # Step 1: Find all metal atoms and categorize residues
    metal_atoms = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                res_id = residue.get_id()[1]
                
                if res_name in standard_residues:
                    standard_residues_found.add(f"{res_name}:{res_id}")
                else:
                    nonstandard_residues_found.add(f"{res_name}:{res_id}")
                
                for atom in residue:
                    atom_name = atom.get_name()
                    
                    # Check if atom is a metal
                    is_metal = False
                    for metal in metals:
                        if atom_name == metal or (atom_name.startswith(metal) and atom_name[len(metal):].isdigit()):
                            is_metal = True
                            break
                    
                    if is_metal:
                        metal_atoms.append({
                            'atom': atom,
                            'name': atom_name,
                            'residue': residue,
                            'res_name': res_name,
                            'res_id': res_id,
                            'position': atom.coord
                        })
                        metal_atoms_found.append(f"{atom_name} in {res_name}:{res_id}")
    
    # Step 2: Find standard residues coordinated to metals
    for metal_info in metal_atoms:
        metal_atom = metal_info['atom']
        metal_name = metal_info['name']
        metal_res_id = metal_info['res_id']
        metal_position = metal_info['position']
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in standard_residues:
                        res_id = residue.get_id()[1]
                        
                        # Skip if this is the same residue as the metal (should not happen with standard residues)
                        if res_id == metal_res_id:
                            continue
                        
                        for atom in residue:
                            distance = np.linalg.norm(metal_position - atom.coord)
                            
                            if distance <= distance_cutoff:
                                
                                metal_coordinated_residues.add(res_id)
                                standard_residues_to_process.add(res_id)
                                
                                coordination_data.append(
                                    f"bond PRO.{metal_res_id}.{metal_name} "
                                    f"PRO.{res_id}.{atom.get_name()}"
                                )
    
    # Step 3: Find non-standard residues containing metals
    metal_containing_nonstandard = {}  # Maps residue ID to the metal info
    
    for metal_info in metal_atoms:
        res_name = metal_info['res_name']
        res_id = metal_info['res_id']
        
        # If this is a non-standard residue
        if res_name not in standard_residues:
            if res_id not in metal_containing_nonstandard:
                metal_containing_nonstandard[res_id] = []
            metal_containing_nonstandard[res_id].append(metal_info)
    
    # Step 4: Find standard residues that bond with non-standard residues containing metals
    # Collect all heavy atoms from standard residues and non-standard residues with metals
    standard_heavy_atoms = []
    standard_atom_info = {}  # Maps atom to (res_id, atom_name)
    
    nonstandard_heavy_atoms = []
    nonstandard_atom_info = {}  # Maps atom to (res_id, atom_name)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()
                
                # If this is a standard residue
                if res_name in standard_residues:
                    for atom in residue:
                        element = atom.element if hasattr(atom, 'element') else atom.get_name()[0]
                        if element != 'H':  # Only heavy atoms
                            standard_heavy_atoms.append(atom)
                            standard_atom_info[atom] = (res_id, atom.get_name())
                
                # If this is a non-standard residue with a metal
                elif res_id in metal_containing_nonstandard:
                    for atom in residue:
                        # Skip the metal atoms themselves
                        is_metal_atom = False
                        for metal_info in metal_containing_nonstandard[res_id]:
                            if atom == metal_info['atom']:
                                is_metal_atom = True
                                break
                        
                        if not is_metal_atom:  # Only include non-metal atoms
                            element = atom.element if hasattr(atom, 'element') else atom.get_name()[0]
                            if element != 'H':  # Only heavy atoms
                                nonstandard_heavy_atoms.append(atom)
                                nonstandard_atom_info[atom] = (res_id, atom.get_name())
     
    # Find bonds between non-standard residues with metals and standard residues
    if nonstandard_heavy_atoms and standard_heavy_atoms:
        # Create neighbor search for standard residue atoms
        ns = NeighborSearch(standard_heavy_atoms)
        
        nonstandard_residues_linked_to_standard = set()
        
        for atom in nonstandard_heavy_atoms:
            nonstandard_res_id, atom_name = nonstandard_atom_info[atom]
            
            # Find nearby standard residue atoms
            nearby_atoms = ns.search(atom.coord, bond_cutoff, level='A')
            
            for nearby_atom in nearby_atoms:
                if nearby_atom in standard_atom_info:
                    standard_res_id, standard_atom_name = standard_atom_info[nearby_atom]
                    
                    
                    nonstandard_standard_bond_data.append(
                        f"bond PRO.{nonstandard_res_id}.{atom_name} "
                        f"PRO.{standard_res_id}.{standard_atom_name}"
                    )
                    
                    # Add this standard residue to the list to process for peptide bonds
                    standard_residues_to_process.add(standard_res_id)
                    nonstandard_residues_linked_to_standard.add(nonstandard_res_id)
    
    # Step 5: Find peptide bonds for all standard residues that are either:
    # 1. Coordinated to a metal atom
    # 2. Bonded to a non-standard residue containing a metal
    # Using a set to track unique bonds and avoid duplicates
    peptide_bonds_set = set()
    
    for model in structure:
        for chain in model:
            chain_residues = list(chain)
            
            for i, residue in enumerate(chain_residues):
                res_id = residue.get_id()[1]
                
                if res_id in standard_residues_to_process:
                    # Check previous residue (N bond)
                    if i > 0:
                        prev_residue = chain_residues[i-1]
                        prev_res_id = prev_residue.get_id()[1]
                        
                        # Create a canonical representation of the bond
                        # Sort the residue IDs to ensure consistent representation
                        bond_parts = sorted([
                            f"PRO.{res_id}.N", 
                            f"PRO.{prev_res_id}.C"
                        ])
                        bond = f"bond {bond_parts[0]} {bond_parts[1]}"
                        peptide_bonds_set.add(bond)
                    
                    # Check next residue (C bond)
                    if i < len(chain_residues) - 1:
                        next_residue = chain_residues[i+1]
                        next_res_id = next_residue.get_id()[1]
                        
                        # Create a canonical representation of the bond
                        # Sort the residue IDs to ensure consistent representation
                        bond_parts = sorted([
                            f"PRO.{res_id}.C", 
                            f"PRO.{next_res_id}.N"
                        ])
                        bond = f"bond {bond_parts[0]} {bond_parts[1]}"
                        peptide_bonds_set.add(bond)
    
    # Convert to list for further processing
    peptide_bond_data = list(peptide_bonds_set)
    
    # Combine and write results to the output file
    all_data = coordination_data + nonstandard_standard_bond_data + peptide_bond_data
    
    if all_data:
        with open(output_file, "w") as f:
            f.write("\n".join(all_data))
    else:
        print("No coordination data found.")

if __name__ == "__main__":
    input_pdb = "easyPARM.pdb"
    
    batch_update_mol2_files('charges_all.dat')
    generate_standard_residue_coordination(input_pdb, output_file="coordinated_residues.txt", distance_cutoff=2.5)
