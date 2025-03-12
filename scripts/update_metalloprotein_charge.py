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
#                              |  $$$$$$/              Ver. 3.10 - 12 February 2025                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


import re
from Bio.PDB import PDBParser
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
                                           distance_cutoff=2.5):
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
    metal_coordinated_residues = set()
    
    # First pass: find metal-coordinated residues
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name().startswith(tuple(metals)):
                        metal_atom = atom
                        metal_name = atom.get_name().split()[0]
                        metal_residue_number = residue.get_id()[1]
                        metal_position = metal_atom.coord
                        
                        # Check for coordinating residues
                        for chain2 in model:
                            for residue2 in chain2:
                                if residue2.get_resname() not in standard_residues:
                                    continue
                                
                                for atom2 in residue2:
                                    distance = np.linalg.norm(metal_position - atom2.coord)
                                    
                                    if distance <= distance_cutoff:
                                        coordinated_residue_number = residue2.get_id()[1]
                                        
                                        # Store coordinated residue numbers
                                        metal_coordinated_residues.add(coordinated_residue_number)
                                        
                                        coordination_data.append(
                                            f"bond PRO.{metal_residue_number}.{metal_name} "
                                            f"PRO.{coordinated_residue_number}.{atom2.get_name()}"
                                        )
    
    # Second pass: find peptide bonds for metal-coordinated residues
    for model in structure:
        for chain in model:
            # Convert chain to list for easier neighbor access
            chain_residues = list(chain)
            
            for i, residue in enumerate(chain_residues):
                current_residue_number = residue.get_id()[1]
                
                # Only process if current residue is metal-coordinated
                if current_residue_number in metal_coordinated_residues:
                    # Check previous residue (N bond)
                    if i > 0:
                        prev_residue = chain_residues[i-1]
                        peptide_bond_data.append(
                            f"bond PRO.{current_residue_number}.N PRO.{prev_residue.get_id()[1]}.C"
                        )
                    
                    # Check next residue (C bond)
                    if i < len(chain_residues) - 1:
                        next_residue = chain_residues[i+1]
                        peptide_bond_data.append(
                            f"bond PRO.{current_residue_number}.C PRO.{next_residue.get_id()[1]}.N"
                        )
    
    # Combine and write results to the output file
    all_data = coordination_data + peptide_bond_data
    
    if all_data:
        with open(output_file, "w") as f:
            f.write("\n".join(all_data))
    else:
        print("No coordination data found.")

if __name__ == "__main__":
    input_pdb = "easyPARM.pdb"
    
    batch_update_mol2_files('charges_all.dat')
    generate_standard_residue_coordination(input_pdb, output_file="coordinated_residues.txt", distance_cutoff=2.5)
