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
import periodictable

# Read the pair atoms and their bond type (single, double, triple)
def read_distance_type_data(distance_type_filename):
    distance_type_data = []
    with open(distance_type_filename, 'r') as distance_type_file:
        for line in distance_type_file:
            parts = line.strip().split()
            distance_type_data.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return distance_type_data

# update the section of connectivity and atomic number in library file
def parse_and_update_complex_lib(filename, distance_type_filename):
    with open(filename, 'r') as file:
        content = file.read()

    distance_type_data = read_distance_type_data(distance_type_filename)

    # Extract and update atoms section
    atoms_pattern = r'(!entry\.mol\.unit\.atoms table.*?!entry\.mol\.unit\.atomspertinfo)'
    atoms_match = re.search(atoms_pattern, content, re.DOTALL)
    
    if atoms_match:
        original_atoms_section = atoms_match.group(1)
        updated_atoms_section = update_atoms_section(original_atoms_section)
        content = content.replace(original_atoms_section, updated_atoms_section)
    else:
        print("Atoms section not found in the file.")

    # Extract and update connectivity section
    connectivity_pattern = r'(!entry\.mol\.unit\.connectivity.*?!entry\.mol\.unit\.hierarchy)'
    connectivity_match = re.search(connectivity_pattern, content, re.DOTALL)
    
    if connectivity_match:
        original_connectivity_section = connectivity_match.group(1)
        updated_connectivity_section = update_connectivity_section(original_connectivity_section, distance_type_data)
        content = content.replace(original_connectivity_section, updated_connectivity_section)
    else:
        print("Connectivity section not found in the file.")

    # Write the updated content back to the file
    with open(filename, 'w') as file:
        file.write(content)

# Keeping the space and read the atom to detect its atomic number
def update_atoms_section(section):
    lines = section.split('\n')
    updated_lines = []
    for line in lines:
        if re.match(r'\s*"', line):
            leading_spaces = re.match(r'(\s*)', line).group(1)
            parts = line.split()
            if len(parts) >= 8:
                atom_label = parts[0].strip('"')
                element = re.sub(r'\d+$', '', atom_label)
                
                try:
                    correct_atomic_number = int(periodictable.elements.symbol(element).number)
                except ValueError:
                    print(f"Warning: Unknown element {element}")
                    updated_lines.append(line)
                    continue

                given_atomic_number = int(parts[6])

                if given_atomic_number != correct_atomic_number:
                    parts[6] = str(correct_atomic_number)
                    updated_line = leading_spaces + ' '.join(parts)
                    updated_lines.append(updated_line)
                else:
                    updated_lines.append(line)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    return '\n'.join(updated_lines)

# update_connectivity_section to only add new bonds
def update_connectivity_section(section, distance_type_data):
    lines = section.split('\n')
    updated_lines = []
    existing_connections = set()
    last_connectivity_line_index = -1

    # First pass: preserve all existing connections without modification
    for i, line in enumerate(lines):
        updated_lines.append(line)  # Keep the line exactly as is
        if line.strip().startswith('!entry'):
            continue
        elif re.match(r'\s*\d', line):
            parts = line.split()
            if len(parts) >= 3:
                atom1, atom2 = map(int, parts[:2])
                # Add to existing connections set (always store smaller index first)
                existing_connections.add((min(atom1, atom2), max(atom1, atom2)))
                last_connectivity_line_index = i

    # Second pass: add only new connections from distance_type_data
    if last_connectivity_line_index != -1:
        # Get the leading spaces from the last connectivity line for consistent formatting
        leading_spaces = re.match(r'(\s*)', lines[last_connectivity_line_index]).group(1)
        
        # Add only bonds that don't exist yet
        for atom1, atom2, bond_type in distance_type_data:
            if (min(atom1, atom2), max(atom1, atom2)) not in existing_connections:
                new_line = f"{leading_spaces}{atom1} {atom2} {bond_type}"
                updated_lines.insert(last_connectivity_line_index + 1, new_line)
                last_connectivity_line_index += 1

    return '\n'.join(updated_lines)

# Input
lib_filename = 'COMPLEX.lib'
distance_type_filename = 'distance_type.dat'
parse_and_update_complex_lib(lib_filename, distance_type_filename)

