import os
import re


#Check if the atom name corresponds to a metal or an element with a large atomic number
# Extract the base element symbol (remove numbers)
def is_metal_or_large_atomic_number(atom_name):
    base_symbol = re.match(r'([A-Za-z]+)', atom_name).group(1).lower()

    periodic_table = {
        'al': 13, 'si': 14, 'ar': 18, 'k': 19, 'cr': 24,  
        'sc': 21, 'ti': 22, 'v': 23, 'mn': 25, 'fe': 26, 'co': 27, 'ni': 28, 
        'cu': 29, 'zn': 30, 'ga': 31, 'ge': 32, 'as': 33, 'se': 34, 'kr': 36, 'rb': 37,
        'sr': 38, 'y': 39, 'zr': 40, 'nb': 41, 'mo': 42, 'tc': 43, 'ru': 44, 'rh': 45, 'pd': 46,
        'ag': 47, 'cd': 48, 'in': 49, 'sn': 50, 'sb': 51, 'te': 52, 'i': 53, 'xe': 54,
        'cs': 55, 'ba': 56, 'la': 57, 'ce': 58, 'pr': 59, 'nd': 60, 'pm': 61, 'sm': 62, 'eu': 63, 'gd': 64,
        'tb': 65, 'dy': 66, 'ho': 67, 'er': 68, 'tm': 69, 'yb': 70, 'lu': 71, 'hf': 72, 'ta': 73,
        'w': 74, 're': 75, 'os': 76, 'ir': 77, 'pt': 78, 'au': 79, 'hg': 80, 'tl': 81, 'pb': 82,
        'bi': 83, 'po': 84, 'at': 85, 'rn': 86, 'fr': 87, 'ra': 88, 'ac': 89, 'th': 90, 'pa': 91,
        'u': 92, 'np': 93, 'pu': 94, 'am': 95, 'cm': 96, 'bk': 97, 'cf': 98, 'es': 99, 'fm': 100, 'md': 101, 'no': 102, 'lr': 103
    }
    # Get the atomic number, defaulting to 0 if not found
    atomic_number = periodic_table.get(base_symbol, 0)
    # Return True if atomic number is 5 or greater
    return atomic_number >= 5 if atomic_number > 0 else False

# Replace metal atoms in the mol2 file with placeholder atoms
def replace_metals(mol2_file, output_file, replacements):
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

            if atom_section:
                parts = line.split()
                if len(parts) > 5:
                    atom_name = parts[1]  # Now checking the atom name
                    if is_metal_or_large_atomic_number(atom_name):
                        # Replace the metal atom with a placeholder
                        replacement = replacements[metal_count % len(replacements)]
                        replacement_map[replacement] = parts[5]  # Store the replacement based on atom type
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

    # Check if more than 8 atoms were replaced
    if len(replaced_atoms) > 6:
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
    replacements = ['s6','p4', 'ss', 'p3','sh', 'p5', 's4', 'p2', 's2', 'SH', 'P', 'S1' ,'N3', 'ce', 'cq', 'c6', 'cz', 'cx','no'   ]

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
