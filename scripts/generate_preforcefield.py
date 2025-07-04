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
#                              |  $$$$$$/              Ver. 4.00 - 8 June 2025                                    #
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


#Returns the periodic table dictionary with or without boron based on the parameter.
#By default includes boron unless metals_complete.dat exists.
def get_periodic_table(include_boron=True):  # Changed default to True
    periodic_table = {
        'al': 13, 'ar': 18, 'cr': 24,
        'sc': 21, 'ti': 22, 'v': 23, 'mn': 25, 'fe': 26,
        'co': 27, 'ni': 28, 'cu': 29, 'zn': 30, 'ga': 31,
        'se': 34, 'kr': 36, 'rb': 37,
        'sr': 38, 'y': 39, 'zr': 40, 'nb': 41, 'mo': 42,
        'tc': 43, 'ru': 44, 'rh': 45, 'pd': 46, 'ag': 47,
        'cd': 48, 'in': 49, 'sn': 50, 'xe': 54, 'cs': 55, 'ba': 56, 'la': 57,
        'ce': 58, 'pr': 59, 'nd': 60, 'pm': 61, 'sm': 62,
        'eu': 63, 'gd': 64, 'tb': 65, 'dy': 66, 'ho': 67,
        'er': 68, 'tm': 69, 'yb': 70, 'lu': 71, 'hf': 72,
        'ta': 73, 'w': 74, 're': 75, 'os': 76, 'ir': 77,
        'pt': 78, 'au': 79, 'hg': 80, 'tl': 81, 'pb': 82,
        'bi': 83, 'po': 84, 'at': 85, 'rn': 86, 'fr': 87,
        'ra': 88, 'ac': 89, 'th': 90, 'pa': 91, 'u': 92,
        'np': 93, 'pu': 94, 'am': 95, 'cm': 96, 'bk': 97,
        'cf': 98, 'es': 99, 'fm': 100, 'md': 101, 'no': 102,
        'lr': 103
    }
    
    if include_boron:
        periodic_table['b'] = 5
        periodic_table['si'] = 14
        periodic_table['ge'] = 32
        periodic_table['as'] = 33
        periodic_table['sb'] = 51
        periodic_table['te'] = 52
        
    return periodic_table

#Determines if an atom is a metal or has a large atomic number.
#Includes boron by default unless metals_complete.dat exists.
def is_metal_or_large_atomic_number(atom_name):
    # Check if metals_complete.dat exists - if it exists, don't include boron
    include_boron = not os.path.exists('metals_complete.dat')
    
    # Get the base symbol from the atom name
    base_symbol = re.match(r'([A-Za-z]+)', atom_name).group(1).lower()
    
    # Get the appropriate periodic table based on file existence
    periodic_table = get_periodic_table(include_boron)
    
    # Get the atomic number, defaulting to 0 if not found
    atomic_number = periodic_table.get(base_symbol, 0)
    
    # Return True if atomic number is 5 or greater
    return atomic_number >= 5 if atomic_number > 0 else False

#Extract all existing atom types from the mol2 file.
def get_existing_atom_types(mol2_file):
    existing_types = set()
    
    try:
        with open(mol2_file, 'r') as file:
            lines = file.readlines()
        
        atom_section = False
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                atom_section = True
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                atom_section = False
                break
            
            if atom_section and line.strip():
                parts = line.split()
                if len(parts) > 5:
                    atom_type = parts[5].strip()
                    existing_types.add(atom_type)
    
    except FileNotFoundError:
        print(f"Warning: {mol2_file} not found. Proceeding with empty existing types set.")
    
    return existing_types


# Comprehensive list of potential replacement atom types
def get_available_replacements(existing_atom_types):
    all_replacements = [
        's6', 'p3', 'ss', 'p5', 's4', 'p2', 'sh', 'S1', 'N3', 
        'ce', 'cq', 'c6', 'cz', 'cx', 'no', 'n1', 'n2', 'n4', 'na', 'nb', 
        'nc', 'nd', 'ne', 'nf', 'nh', 'ni', 'nj', 'nk', 'nl', 'nm', 'nn', 
        'np', 'nq', 'nr', 'ns', 'nt', 'nu', 'nv', 'nw', 'nx', 'ny', 'nz',
        'c1', 'c2', 'c3', 'c4', 'c5', 'c7', 'c8', 'c9', 'ca', 'cb', 'cc', 
        'cd', 'cf', 'cg', 'ch', 'ci', 'cj', 'ck', 'cl', 'cm', 'cn', 'co', 
        'cp', 'cr', 'cs', 'ct', 'cu', 'cv', 'cw', 'cy', 'cA', 'cB', 'cC',
        'o1', 'o2', 'o3', 'o4', 'o5', 'o6', 'o7', 'o8', 'o9', 'oa', 'ob', 
        'oc', 'od', 'oe', 'of', 'og', 'oh', 'oi', 'oj', 'ok', 'ol', 'om', 
        'on', 'oo', 'op', 'oq', 'or', 'os', 'ot', 'ou', 'ov', 'ow', 'ox',
        's1', 's3', 's5', 's7', 's8', 's9', 'sa', 'sb', 'sc', 'sd', 'se', 
        'sf', 'sg', 'si', 'sj', 'sk', 'sl', 'sm', 'sn', 'so', 'sp', 'sq', 
        'sr', 'st', 'su', 'sv', 'sw', 'sx', 'sy', 'sz'
    ]
    
    # Filter out atom types that already exist in the mol2 file
    available_replacements = [atom_type for atom_type in all_replacements 
                             if atom_type not in existing_atom_types]
    if not available_replacements:
        raise ValueError("No available replacement atom types found! All potential replacements already exist in the mol2 file.")
    
    return available_replacements


#Replace metal atoms in the mol2 file with placeholder atoms that don't conflict with existing types.
def replace_metals(mol2_file, output_file, existing_atom_types=None):
    # Get existing atom types if not provided
    if existing_atom_types is None:
        existing_atom_types = get_existing_atom_types(mol2_file)
    
    # Get available replacement atom types
    available_replacements = get_available_replacements(existing_atom_types)
    
    replacement_map = {}  # To keep track of what was replaced by what
    replaced_atoms = []  # To store the atom numbers of replaced atoms

    with open(mol2_file, 'r') as file:
        lines = file.readlines()

    atom_section = False
    metal_count = 0
    
    with open(output_file, 'w') as outfile:
        for line in lines:
            # Identify the atom section
            if line.startswith('@<TRIPOS>ATOM'):
                atom_section = True
                outfile.write(line)
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                atom_section = False

            if atom_section and line.strip():
                parts = line.split()
                if len(parts) > 5:
                    atom_name = parts[1]  # Now checking the atom name
                    if is_metal_or_large_atomic_number(atom_name):
                        # Check if we have enough available replacements
                        if metal_count >= len(available_replacements):
                            raise ValueError(f"Not enough available replacement atom types! "
                                           f"Found {metal_count + 1} metals but only {len(available_replacements)} available replacements.")
                        
                        # Replace the metal atom with an available placeholder
                        replacement = available_replacements[metal_count]
                        original_atom_type = parts[5]
                        replacement_map[replacement] = original_atom_type  # Store the replacement mapping
                        replaced_atoms.append(parts[0])  # Save the atom number
                        parts[5] = replacement
                        metal_count += 1
                        
                        # Format the line with the new replacement
                        line = f"{parts[0]:>7} {parts[1]:<8} {parts[2]:>10} {parts[3]:>10} {parts[4]:>10} {parts[5]:<8} {parts[6]:>3} {parts[7]:<6} {parts[8]:>10}\n"
            
            outfile.write(line)

    # Save the replaced atom numbers to a file
    with open('metal_number.dat', 'w') as metal_file:
        for atom_number in replaced_atoms:
            metal_file.write(f"{atom_number}\n")

    # Check if more than 28 atoms were replaced
    if len(replaced_atoms) > 28:
        with open('more_metal.dat', 'w') as more_metal_file:
            for atom_number in replaced_atoms:
                more_metal_file.write(f"{atom_number}\n")

    return replacement_map

#Reverse the metal replacements in the mol2 file.
def reverse_replace_mol2(mol2_file, replacement_map):
    atom_types = {}  # To store the original atom types

    with open(mol2_file, 'r') as file:
        lines = file.readlines()

    atom_section = False
    with open(mol2_file, 'w') as outfile:
        for line in lines:
            if line.startswith('@<TRIPOS>ATOM'):
                atom_section = True
                outfile.write(line)
                continue
            elif line.startswith('@<TRIPOS>BOND'):
                atom_section = False

            if atom_section:
                parts = line.split()
                if len(parts) > 5:
                    atom_symbol = parts[5]
                    if atom_symbol in replacement_map:
                        # Reverse the replacement
                        original_atom = replacement_map[atom_symbol]
                        parts[5] = original_atom
                        atom_types[parts[0]] = original_atom  # Store atom number and type
                        line = f"{parts[0]:>7} {parts[1]:<8} {parts[2]:>10} {parts[3]:>10} {parts[4]:>10} {parts[5]:<8} {parts[6]:>3} {parts[7]:<6} {parts[8]:>10}\n"
            outfile.write(line)

    # Save the atom types to atom_type.dat
    with open('atom_type.dat', 'w') as atom_type_file:
        for atom_number, atom_type in atom_types.items():
            atom_type_file.write(f"{atom_type}\n")

#Reverse the metal replacements in the frcmod file
def reverse_replace_frcmod(frcmod_file, replacement_map):
    # Define the section headers
    section_headers = ['MASS', 'BOND', 'ANGLE', 'DIHE', 'IMPROPER', 'NONBON']
    
    with open(frcmod_file, 'r') as file:
        lines = file.readlines()

    current_section = None
    with open(frcmod_file, 'w') as outfile:
        for line in lines:
            # Check if the line starts a new section
            if any(line.startswith(header) for header in section_headers):
                current_section = line.strip()
                outfile.write(line)
                continue

            # Only perform replacements if the current section is between the defined headers
            if current_section in section_headers:
                original_line = line  # Keep the original line for reference
                for replacement, original_symbol in replacement_map.items():
                    if replacement in line:
                        # Replace the element ensuring it takes up exactly 4 characters
                        index = line.find(replacement)
                        if index != -1:
                            original_symbol = original_symbol.ljust(2)  # Ensure the original symbol is padded to 4 characters
                            line = line[:index] + original_symbol + line[index + len(replacement):]
            outfile.write(line)

def main():
    input_file = 'COMPLEX.mol2'
    output_file = 'COMPLEX_modified.mol2'
    frcmod_file = 'COMPLEX_modified.frcmod'
    
    # Define the list of replacement atoms
    replacements = get_existing_atom_types(input_file) 

    # Replace metals in the mol2 file and store the replacements
    replacement_map = replace_metals(input_file, output_file, replacements)

    # Check if the .frcmod file exists
    if not os.path.exists(frcmod_file):
        return

    # Reverse replacements in the mol2 file
    reverse_replace_mol2(output_file, replacement_map)
    # Reverse replacements in the frcmod file
    reverse_replace_frcmod(frcmod_file, replacement_map)

if __name__ == "__main__":
    main()

