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


# Function to extract bond data from a mol2 file
def extract_bond_data(mol2_filename):
    bond_data = []
    sequence_number = 0
    with open(mol2_filename, 'r') as mol2_file:
        is_bond_section = False
        for line in mol2_file:
            line = line.strip()
            # Identify the bond section
            if line == "@<TRIPOS>BOND":
                is_bond_section = True
                continue
            if line == "@<TRIPOS>SUBSTRUCTURE":
                is_bond_section = False
            # Process lines in the bond section
            if is_bond_section and line:
                parts = line.split()
                try:
                    # Handle both integer and 'ar' cases for parts[3]
                    if parts[3] == 'ar':
                        bond_data.append((int(parts[1]), int(parts[2]), 'ar'))
                    else:
                        bond_data.append((int(parts[1]), int(parts[2]), int(parts[3])))
                    
                    sequence_number = int(parts[0])  # Keep track of the last sequence number
                except ValueError as e:
                    print(f"Error processing line: {line} - {e}")
    return bond_data, sequence_number

# Function to read distance data from a file
def read_distance_type_data(distance_type_filename):
    distance_type_data = []
    with open(distance_type_filename, 'r') as distance_type_file:
        for line in distance_type_file:
            parts = line.strip().split()
            distance_type_data.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return distance_type_data

# Function to get the element symbol from an atom name
def get_element_symbol(atom_name):
    # Strip any digits from the atom name to get the element symbol
    for i, char in enumerate(atom_name):
        if char.isdigit():
            return atom_name[:i].capitalize()
    return atom_name.capitalize()

# Function to replace 'DU' or 'du' with the correct atom symbol in a mol2 file
def replace_du_with_atom_symbol(mol2_filename):
    with open(mol2_filename, 'r') as mol2_file:
        lines = mol2_file.readlines()
    
    # Find the atom and bond sections
    atom_section_start = None
    bond_section_start = None
    
    for i, line in enumerate(lines):
        if "@<TRIPOS>ATOM" in line:
            atom_section_start = i
        elif "@<TRIPOS>BOND" in line:
            bond_section_start = i
            break  # No need to keep searching after finding the bond section
    
    if atom_section_start is None or bond_section_start is None:
        raise ValueError("Could not find atom or bond section in the file.")
    
    # Modify the lines in the ATOM section
    for i in range(atom_section_start + 1, bond_section_start):
        parts = lines[i].split()
        if len(parts) > 5 and parts[5].lower() == "du":
            atom_name = parts[1]
            element_symbol = get_element_symbol(atom_name)
            # Replace 'DU' or 'du' with the correct atom symbol, preserving the format
            start_idx = lines[i].find(parts[5])
            lines[i] = lines[i][:start_idx] + element_symbol.title() + lines[i][start_idx + len(parts[5]):]
    
    # Write the modified lines back to the file
    with open(mol2_filename, 'w') as mol2_file:
        mol2_file.writelines(lines)

# Function to replace the bond count in the mol2 file header
def replace_bond_count_in_mol2(mol2_filename, total_bonds):
    with open(mol2_filename, 'r') as mol2_file:
        lines = mol2_file.readlines()
    
    if len(lines) < 3:
        raise ValueError("The mol2 file is too short to contain the required header lines.")

    # Modify the third line to update the bond count, preserving the original format
    header_line = lines[2].strip()
    header_parts = header_line.split()
    
    if len(header_parts) >= 5:  # Ensure we have at least 5 spaces
        # Replace the second number (bond count) with the new total_bonds
        header_parts[1] = str(total_bonds)
        
        # Reconstruct the line, preserving original spacing
        new_header_line = ''
        current_index = 0
        for part in header_parts:
            # Find the next non-space character
            while current_index < len(header_line) and header_line[current_index].isspace():
                new_header_line += header_line[current_index]
                current_index += 1
            
            # Add the new part
            new_header_line += part
            
            # Move the index past this part
            current_index += len(header_line[current_index:].split(None, 1)[0])
        
        # Add any trailing whitespace
        new_header_line += header_line[current_index:]
        
        # Ensure the line ends with a newline
        if not new_header_line.endswith('\n'):
            new_header_line += '\n'
        
        lines[2] = new_header_line
    else:
        raise ValueError("Unexpected format in the third line of the mol2 file.")

    # Write the modified lines back to the file
    with open(mol2_filename, 'w') as mol2_file:
        mol2_file.writelines(lines)

def update_bonds_in_mol2(mol2_filename, updated_bonds, new_bonds):
    with open(mol2_filename, 'r') as mol2_file:
        lines = mol2_file.readlines()

    # Find the bond and substructure sections
    bond_section_start = None
    bond_section_end = None

    for i, line in enumerate(lines):
        if "@<TRIPOS>BOND" in line:
            bond_section_start = i
        elif "@<TRIPOS>SUBSTRUCTURE" in line:
            bond_section_end = i
            break

    if bond_section_start is None or bond_section_end is None:
        raise ValueError("Could not find bond or substructure section in the file.")

    # Extract existing bonds without modification
    existing_bonds = [line for line in lines[bond_section_start + 1:bond_section_end] if line.strip()]

    # Add new bonds
    last_bond_number = int(existing_bonds[-1].split()[0]) if existing_bonds else 0
    for atom1, atom2, bond_type in new_bonds:
        last_bond_number += 1
        existing_bonds.append(f"{last_bond_number:>5} {atom1:>5} {atom2:>5} {bond_type}\n")

    # Write the modified file with existing bonds (unchanged) and new bonds
    with open(mol2_filename, 'w') as mol2_file:
        mol2_file.writelines(lines[:bond_section_start + 1])  # Write up to @<TRIPOS>BOND
        mol2_file.writelines(existing_bonds)  # Write all bonds (existing and new)
        mol2_file.writelines(lines[bond_section_end:])  # Write from @<TRIPOS>SUBSTRUCTURE to end

def check_and_extend_bond_data(mol2_filename, distance_type_filename):
    # Extract existing bond data and last sequence number
    bond_data, _ = extract_bond_data(mol2_filename)
    # Read distance and bond type data
    distance_type_data = read_distance_type_data(distance_type_filename)

    # Compare and update bonds
    bond_dict = {(atom1, atom2): bond_type for atom1, atom2, bond_type in bond_data}
    updated_bonds = {}
    new_bonds = []
    
    for atom1, atom2, new_bond_type in distance_type_data:
        if (atom1, atom2) in bond_dict:
            if bond_dict[(atom1, atom2)] != new_bond_type:
                updated_bonds[(atom1, atom2)] = new_bond_type
        elif (atom2, atom1) in bond_dict:
            if bond_dict[(atom2, atom1)] != new_bond_type:
                updated_bonds[(atom2, atom1)] = new_bond_type
        else:
            new_bonds.append((atom1, atom2, new_bond_type))
    
    total_bonds = len(bond_data) + len(new_bonds)
    
    # Replace DU/du with the correct atom symbol in the mol2 file
    replace_du_with_atom_symbol(mol2_filename)
    
    # Replace the bond count in the mol2 file's third line
    replace_bond_count_in_mol2(mol2_filename, total_bonds)
    
    # Update existing bonds and add new ones
    update_bonds_in_mol2(mol2_filename, updated_bonds, new_bonds)

mol2_filename = "COMPLEX.mol2"
distance_type_filename = "distance_type.dat"
check_and_extend_bond_data(mol2_filename, distance_type_filename)# Function to insert missing bonds into the mol2 file
