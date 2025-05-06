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
#                              |  $$$$$$/              Ver. 3.30 - 5 May 2025                                     #
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

# Dictionary of atomic numbers for elements
ATOMIC_NUMBERS = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
    'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
    'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
    'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
    'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
    'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
    'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
    'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
    'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
    'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
    'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
    'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}

def get_atomic_number(atom_name):
    element = ''.join(c for c in atom_name if c.isalpha())
    element = element.capitalize()
    return ATOMIC_NUMBERS.get(element, 0)

def extract_atom_info_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_info = {}
        atom_types_present = set()  # New set to track all atom types
        is_atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    atom_type = parts[5].split('.')[0] if '.' in parts[5] else parts[5]  # Extract atom type
                    atom_types_present.add(atom_type)  # Add to set of present types
                    atom_name_base = re.sub(r'\d+', '', atom_name)
                    atomic_num = get_atomic_number(atom_name_base)
                    if atomic_num > 10:
                        atom_info[atom_id] = atom_name_base
        return atom_info, atom_types_present  # Return both the atom info and the set of present types
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}, set()

def read_metal_numbers(file_path, atom_info):
    try:
        with open(file_path, 'r') as file:
            numbers = [int(line.strip()) for line in file if line.strip()]
        return [num for num in numbers if num in atom_info]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []

def count_bonds(metal_numbers, distance_file_path):
    bond_counts = {metal: 0 for metal in metal_numbers}
    try:
        with open(distance_file_path, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) >= 2:
                    atom1 = int(parts[0])
                    atom2 = int(parts[1])
                    if atom1 in bond_counts:
                        bond_counts[atom1] += 1
                    if atom2 in bond_counts:
                        bond_counts[atom2] += 1
        return bond_counts
    except Exception as e:
        print(f"Error reading {distance_file_path}: {str(e)}")
        return {}

def read_uff_data():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        uff_data_file_path = os.path.join(script_dir, 'uff_data.txt')
        
        with open(uff_data_file_path, 'r') as file:
            uff_data = {}
            for line in file:
                parts = line.split()
                if len(parts) >= 7:
                    atom_type = parts[0]
                    col4_half = float(parts[3]) / 2.0
                    col5 = float(parts[4])
                    uff_data[atom_type] = (col4_half, col5)
        return uff_data
    except Exception as e:
        print(f"Error reading {uff_data_file_path}: {str(e)}")
        return {}

def update_frcmod_with_uff_data(frcmod_file_path, uff_data, atom_info, bond_counts, atom_types_present):
    try:
        with open(frcmod_file_path, 'r') as file:
            lines = file.readlines()

        nonbon_start_index = None
        for i, line in enumerate(lines):
            if line.strip().startswith("NONBON"):
                nonbon_start_index = i + 1
                break

        if nonbon_start_index is None:
            print("No NONBON section found in the frcmod file.")
            return

        # Create new lines list starting with everything up to NONBON section
        new_lines = lines[:nonbon_start_index]
        
        # Process remaining lines
        for line in lines[nonbon_start_index:]:
            parts = line.split()
            if len(parts) >= 3:
                atom_type = parts[0].strip()
                # Only keep the line if the atom type exists in mol2 file
                if atom_type in atom_types_present:
                    # Process metal atoms as before
                    matching_atom_nums = [num for num, aname in atom_info.items() if aname == atom_type]
                    
                    if matching_atom_nums:
                        atom_num = matching_atom_nums[0]
                        if atom_num in bond_counts:
                            bond_count = bond_counts[atom_num]
                            base_atom_type = atom_type

                            matching_uff_keys = [key for key in uff_data.keys() if key.startswith(base_atom_type)]
                            new_atom_type = None

                            for key in matching_uff_keys:
                                if key[len(base_atom_type):].startswith(str(bond_count)):
                                    new_atom_type = key
                                    break

                            if new_atom_type and new_atom_type in uff_data:
                                col2, col3 = uff_data[new_atom_type]
                            else:
                                matched_key = None
                                for key in uff_data.keys():
                                    if key.startswith(base_atom_type) and key != new_atom_type:
                                        matched_key = key
                                        break
                                
                                if matched_key:
                                    col2, col3 = uff_data[matched_key]
                                else:
                                    new_lines.append(line)
                                    continue

                            new_lines.append(f"{atom_type:<10}{col2:>8.4f} {col3:>8.4f}\n")
                        else:
                            new_lines.append(line)
                    else:
                        new_lines.append(line)
            else:
                new_lines.append(line)

        with open(frcmod_file_path, 'w') as file:
            file.writelines(new_lines)

    except Exception as e:
        print(f"Error processing {frcmod_file_path}: {str(e)}")

# Inputs
mol2_file_path = 'NEW_COMPLEX.mol2'
metal_numbers_file_path = 'metal_number.dat'
distance_file_path = 'distance.dat'
frcmod_file_path = 'updated_updated_COMPLEX_modified2.frcmod'

# Run the extraction of info and total bond number
atom_info, atom_types_present = extract_atom_info_from_mol2(mol2_file_path)  # Now gets both returns
metal_numbers = read_metal_numbers(metal_numbers_file_path, atom_info)
bond_counts = count_bonds(metal_numbers, distance_file_path)

# Update the file
uff_data = read_uff_data()
update_frcmod_with_uff_data(frcmod_file_path, uff_data, atom_info, bond_counts, atom_types_present)
