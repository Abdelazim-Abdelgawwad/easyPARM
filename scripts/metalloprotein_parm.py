#!/usr/bin/env python3
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
#                              |  $$$$$$/              Ver. 4.15 - 17 October 2025                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import re
import os
import sys


#Read metalloprotein_atomtype.dat file and return mapping from new_atom_type to original_atom_type
def read_atom_type_mapping(filename):
    new_to_original = {}
    original_to_new = {}
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    new_atom_type = parts[0]
                    original_atom_type = parts[1]
                    new_to_original[new_atom_type] = original_atom_type
                    original_to_new[original_atom_type] = new_atom_type
    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
        exit(1)
        
    return new_to_original, original_to_new


#Read protein_parm.dat file and extract bond and angle parameters from 
def read_protein_parm_data(filename):
    bond_params = {}
    angle_params = {}
    
    try:
        with open(filename, 'r') as f:
            section = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith('BOND'):
                    section = 'BOND'
                    continue
                elif line.startswith('ANGLE'):
                    section = 'ANGLE'
                    continue
                elif line.startswith('DIHE'):
                    section = 'DIHE'
                    continue
                
                if not line or line.startswith('#'):
                    continue
                    
                if section == 'BOND':
                    # Extract bond parameters
                    match = re.match(r'(\S+)-(\S+)\s+(\S+)\s+(\S+)', line)
                    if match:
                        atom1, atom2, force, dist = match.groups()
                        # Store parameters for both directions
                        key1 = f"{atom1}-{atom2}"
                        key2 = f"{atom2}-{atom1}"
                        bond_params[key1] = (float(force), float(dist))
                        bond_params[key2] = (float(force), float(dist))
                        
                elif section == 'ANGLE':
                    # Extract angle parameters
                    match = re.match(r'(\S+)-(\S+)-(\S+)\s+(\S+)\s+(\S+)', line)
                    if match:
                        atom1, atom2, atom3, force, angle = match.groups()
                        # Store parameters for both directions
                        key1 = f"{atom1}-{atom2}-{atom3}"
                        key2 = f"{atom3}-{atom2}-{atom1}"
                        angle_params[key1] = (float(force), float(angle))
                        angle_params[key2] = (float(force), float(angle))
    
    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
        exit(1)
        
    return bond_params, angle_params


#Update COMPLEX.frcmod file with parameters from protein_parm.dat based on atom type mappings.
def update_complex_frcmod(frcmod_filename, new_to_original, original_to_new, protein_bond_params, protein_angle_params):
    try:
        with open(frcmod_filename, 'r') as f:
            frcmod_content = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find file {frcmod_filename}")
        exit(1)
    
    # Process the file content
    section = None
    updated_content = []
    
    for line_num, line in enumerate(frcmod_content):
        original_line = line
        stripped_line = line.strip()
        
        # Determine section
        if stripped_line == 'BOND':
            section = 'BOND'
            updated_content.append(original_line)
            continue
        elif stripped_line == 'ANGLE':
            section = 'ANGLE'
            updated_content.append(original_line)
            continue
        elif stripped_line == 'DIHE':
            section = 'DIHE'
            updated_content.append(original_line)
            continue
        elif not stripped_line:
            updated_content.append(original_line)
            continue
        
        # Process based on section
        if section == 'BOND':
            # Match bond parameters: atom1-atom2 force distance
            match = re.match(r'(\S+)-(\S+)\s+(\S+)\s+(\S+)', stripped_line)
            if match:
                atom1, atom2, force, dist = match.groups()
                
                # Convert new atom types to original atom types
                orig_atom1 = new_to_original.get(atom1, atom1)
                orig_atom2 = new_to_original.get(atom2, atom2)
                
                # Create keys for lookup
                protein_key1 = f"{orig_atom1}-{orig_atom2}"
                protein_key2 = f"{orig_atom2}-{orig_atom1}"
                
                # Check if this bond exists in protein_parm.dat
                if protein_key1 in protein_bond_params:
                    new_force, new_dist = protein_bond_params[protein_key1]
                    # Preserve the original format and spacing
                    prefix_match = re.match(r'(\S+\s+)', line)
                    prefix = prefix_match.group(1) if prefix_match else ""
                    new_line = f"{prefix}{new_force:.1f} {new_dist:.3f}\n"
                    updated_content.append(new_line)
                elif protein_key2 in protein_bond_params:
                    new_force, new_dist = protein_bond_params[protein_key2]
                    # Preserve the original format and spacing
                    prefix_match = re.match(r'(\S+\s+)', line)
                    prefix = prefix_match.group(1) if prefix_match else ""
                    new_line = f"{prefix}{new_force:.1f} {new_dist:.3f}\n"
                    updated_content.append(new_line)
                else:
                    # Keep original line if no match found
                    updated_content.append(original_line)
            else:
                updated_content.append(original_line)
        
        elif section == 'ANGLE':
            # Match angle parameters: atom1-atom2-atom3 force angle
            match = re.match(r'(\S+)-(\S+)-(\S+)\s+(\S+)\s+(\S+)', stripped_line)
            if match:
                atom1, atom2, atom3, force, angle = match.groups()
                
                # Convert new atom types to original atom types
                orig_atom1 = new_to_original.get(atom1, atom1)
                orig_atom2 = new_to_original.get(atom2, atom2)
                orig_atom3 = new_to_original.get(atom3, atom3)
                
                # Create keys for lookup
                protein_key1 = f"{orig_atom1}-{orig_atom2}-{orig_atom3}"
                protein_key2 = f"{orig_atom3}-{orig_atom2}-{orig_atom1}"
                
                # Check if this angle exists in protein_parm.dat
                if protein_key1 in protein_angle_params:
                    new_force, new_angle = protein_angle_params[protein_key1]
                    # Preserve the original format and spacing
                    prefix_match = re.match(r'(\S+\s+)', line)
                    prefix = prefix_match.group(1) if prefix_match else ""
                    new_line = f"{prefix}{new_force:.3f} {new_angle:.3f}\n"
                    updated_content.append(new_line)
                elif protein_key2 in protein_angle_params:
                    new_force, new_angle = protein_angle_params[protein_key2]
                    # Preserve the original format and spacing
                    prefix_match = re.match(r'(\S+\s+)', line)
                    prefix = prefix_match.group(1) if prefix_match else ""
                    new_line = f"{prefix}{new_force:.3f} {new_angle:.3f}\n"
                    updated_content.append(new_line)
                else:
                    # Keep original line if no match found
                    updated_content.append(original_line)
            else:
                updated_content.append(original_line)
        
        else:
            # For sections other than BOND and ANGLE
            updated_content.append(original_line)
    
    # Write updated content back to file
    backup_filename = f"{frcmod_filename}.bak"
    try:
        # Create backup of original file
        with open(backup_filename, 'w') as f:
            for line in frcmod_content:
                f.write(line)
                
        # Write updated content
        with open(frcmod_filename, 'w') as f:
            for line in updated_content:
                f.write(line)
                
    except Exception as e:
        print(f"Error updating file: {e}")
        exit(1)


def main():
    # Check command line arguments
    if len(sys.argv) != 2:
        print("Usage: python3 code.py <parm_file>")
        exit(1)
    
    parm_type = sys.argv[1]
    # Determine which parameter file to use based on argument
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parm_file = os.path.join(script_dir, parm_type)

    if not os.path.exists(parm_type):
        print(f"Error: Parm file '{parm_type}' not found")
        sys.exit(1)
     
    atom_type_file = "metalloprotein_atomtype.dat"
    frcmod_file = "COMPLEX.frcmod"

    # Read atom type mapping
    new_to_original, original_to_new = read_atom_type_mapping(atom_type_file)
    
    # Read protein parameter data
    protein_bond_params, protein_angle_params = read_protein_parm_data(parm_file)
    
    # Update complex frcmod file
    update_complex_frcmod(frcmod_file, new_to_original, original_to_new, protein_bond_params, protein_angle_params)
    
    
if __name__ == "__main__":
    main()
