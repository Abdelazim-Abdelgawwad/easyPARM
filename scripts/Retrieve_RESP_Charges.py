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


import os
import re

#Update charges in a MOL2 file based on a charge file containing only charges.
def update_mol2_file(input_mol2, charge_file, output_mol2=None):
    try:
        # output file
        if output_mol2 is None:
            base, ext = os.path.splitext(input_mol2)
            output_mol2 = f"{base}{ext}"
        
        # Read new atom charges
        new_charges = []
        with open(charge_file, 'r') as f:
            for line in f:
                # each line contains a single charge value
                charge = float(line.strip())
                new_charges.append(charge)
        
        # Process MOL2 file
        updated_mol2_lines = []
        with open(input_mol2, 'r') as input_file:
            is_atom_section = False
            charge_index = 0
            
            for line in input_file:
                # Detect and modify the atom section
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    charge_index = 0
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
                        atom_id, atom_name, coords, atom_type, post_type, old_charge = atom_match.groups()
                        
                        # Update charge if we have new charges available
                        if charge_index < len(new_charges):
                            new_charge = new_charges[charge_index]
                            updated_line = f"{atom_id}{atom_name}{coords}{atom_type}{post_type}{new_charge:8.6f}\n"
                            updated_mol2_lines.append(updated_line)
                            charge_index += 1
                        else:
                            # If we run out of new charges, keep the original line
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

if __name__ == "__main__":
    input_mol2 = "COMPLEX.mol2"
    charge_file = "charges.chg"
    update_mol2_file(input_mol2, charge_file)
