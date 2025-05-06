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

import re
import periodictable
import sys

# Dictionary of special atom labels to their correct element types
SPECIAL_ATOM_MAPPING = {
    "CA": "C", "HA": "H", "CB": "C", "HB1": "H", "HB2": "H", "HB3": "H",
    "CG": "C", "HG2": "H", "HG3": "H", "CD": "C", "HD2": "H", "HD3": "H",
    "NE": "N", "HE": "H", "CZ": "C", "NH1": "N", "HH11": "H", "HH12": "H",
    "NH2": "N", "HH21": "H", "HH22": "H", "OD1": "O", "OD2": "O", "ND2": "N",
    "HD21": "H", "HD22": "H", "SG": "S", "HG": "H", "OE1": "O", "OE2": "O",
    "HE2": "H", "NE2": "N", "HE21": "H", "HE22": "H", "HA2": "H", "HA3": "H",
    "ND1": "N", "HD1": "H", "CE1": "C", "HE1": "H", "CD2": "C", "HD23": "H",
    "HB": "H", "CG2": "C", "HG21": "H", "HG22": "H", "HG23": "H", "CG1": "C",
    "HG12": "H", "HG13": "H", "CD1": "C", "HD11": "H", "HD12": "H", "HD13": "H",
    "CE": "C", "HE3": "H", "NZ": "N", "HZ2": "H", "HZ3": "H", "HZ1": "H",
    "SD": "S", "HZ": "H", "CE2": "C", "OG": "O", "OG1": "O", "HG1": "H",
    "NE1": "N", "CZ2": "C", "CH2": "C", "HH2": "H", "CZ3": "C", "CE3": "C",
    "OH": "O", "HH": "H", "HG11": "H"
}

# Read the pair atoms and their bond type (single, double, triple)
def read_distance_type_data(distance_type_filename):
    distance_type_data = []
    with open(distance_type_filename, 'r') as distance_type_file:
        for line in distance_type_file:
            parts = line.strip().split()
            distance_type_data.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return distance_type_data

# Update the section of connectivity and atomic number in the library file
def parse_and_update_complex_lib(filename, molecule_name, distance_type_filename):
    with open(filename, 'r') as file:
        content = file.read()

    distance_type_data = read_distance_type_data(distance_type_filename)

    # Define patterns dynamically based on the molecule name
    atoms_pattern = rf'(!entry\.{molecule_name}\.unit\.atoms table.*?!entry\.{molecule_name}\.unit\.atomspertinfo)'
    connectivity_pattern = rf'(!entry\.{molecule_name}\.unit\.connectivity.*?!entry\.{molecule_name}\.unit\.hierarchy)'

    # Extract and update atoms section
    atoms_match = re.search(atoms_pattern, content, re.DOTALL)
    if atoms_match:
        original_atoms_section = atoms_match.group(1)
        updated_atoms_section = update_atoms_section(original_atoms_section)
        content = content.replace(original_atoms_section, updated_atoms_section)
        print(f"Updated atoms section for '{molecule_name}'")
    else:
        print(f"Atoms section for '{molecule_name}' not found in the file.")

    # Extract and update connectivity section
    connectivity_match = re.search(connectivity_pattern, content, re.DOTALL)
    if connectivity_match:
        original_connectivity_section = connectivity_match.group(1)
        updated_connectivity_section = update_connectivity_section(original_connectivity_section, distance_type_data)
        content = content.replace(original_connectivity_section, updated_connectivity_section)
        print(f"Updated connectivity section for '{molecule_name}'")
    else:
        print(f"Connectivity section for '{molecule_name}' not found in the file.")

    # Write the updated content back to the file
    with open(filename, 'w') as file:
        file.write(content)

# Process all molecules in the file but only update their atom sections
def parse_and_update_all_molecules_atoms(filename):
    with open(filename, 'r') as file:
        content = file.read()
    
    # Find all molecule names in the file
    molecule_pattern = r'!entry\.([^\.]+)\.unit\.atoms table'
    molecule_names = re.findall(molecule_pattern, content)
    
    # Remove duplicates while preserving order
    unique_molecule_names = []
    for name in molecule_names:
        if name not in unique_molecule_names:
            unique_molecule_names.append(name)
    
    if not unique_molecule_names:
        print(f"No molecule entries found in {filename}")
        return
    
    print(f"Found molecule names: {', '.join(unique_molecule_names)}")
    
    # Process each molecule name
    for molecule_name in unique_molecule_names:
        # Define pattern for this molecule
        atoms_pattern = rf'(!entry\.{molecule_name}\.unit\.atoms table.*?!entry\.{molecule_name}\.unit\.atomspertinfo)'
        
        # Extract and update atoms section
        atoms_match = re.search(atoms_pattern, content, re.DOTALL)
        if atoms_match:
            original_atoms_section = atoms_match.group(1)
            updated_atoms_section = update_atoms_section(original_atoms_section)
            content = content.replace(original_atoms_section, updated_atoms_section)
            print(f"Updated atoms section for '{molecule_name}'")
        else:
            print(f"Atoms section for '{molecule_name}' not found in the file.")
    
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
                
                # First check if this is a special atom label
                if atom_label in SPECIAL_ATOM_MAPPING:
                    element = SPECIAL_ATOM_MAPPING[atom_label]
                else:
                    # Otherwise extract element from label (remove numbers)
                    element = re.sub(r'\d+$', '', atom_label)
                
                try:
                    correct_atomic_number = int(periodictable.elements.symbol(element).number)
                except ValueError:
                    print(f"Warning: Unknown element {element} from atom {atom_label}")
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

# Update connectivity section to only add new bonds
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
if __name__ == "__main__":
    lib_filename = 'COMPLEX.lib'
    molecule_name = 'mol'  # Default molecule name
    distance_type_filename = 'distance_type.dat'
    
    # Process command line arguments if provided
    if len(sys.argv) > 1:
        lib_filename = sys.argv[1]  # First argument: Library file
    
    if len(sys.argv) > 2:
        molecule_name = sys.argv[2]  # Second argument: Molecule name
    
    # Process the specific molecule first (update both atoms and connectivity)
    parse_and_update_complex_lib(lib_filename, molecule_name, distance_type_filename)
    
    # Create a temporary file for this operation to avoid overwriting the changes from the first function
    import tempfile
    import shutil
    
    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_filename = temp_file.name
    temp_file.close()
    
    # Copy the content to the temporary file
    shutil.copy2(lib_filename, temp_filename)
    
    # Read all molecules from the temporary file
    with open(temp_filename, 'r') as file:
        content = file.read()
    
    # Find all molecule names in the file
    molecule_pattern = r'!entry\.([^\.]+)\.unit\.atoms table'
    all_molecule_names = re.findall(molecule_pattern, content)
    
    # Remove duplicates while preserving order
    unique_molecule_names = []
    for name in all_molecule_names:
        if name != molecule_name and name not in unique_molecule_names:
            unique_molecule_names.append(name)
     
    if unique_molecule_names:
        for other_molecule in unique_molecule_names:
            # Define pattern for this molecule
            atoms_pattern = rf'(!entry\.{other_molecule}\.unit\.atoms table.*?!entry\.{other_molecule}\.unit\.atomspertinfo)'
            
            # Extract and update atoms section
            atoms_match = re.search(atoms_pattern, content, re.DOTALL)
            if atoms_match:
                original_atoms_section = atoms_match.group(1)
                updated_atoms_section = update_atoms_section(original_atoms_section)
                content = content.replace(original_atoms_section, updated_atoms_section)
            else:
                print(f"Atoms section for '{other_molecule}' not found in the file.")
        
        # Write the updated content back to the original file
        with open(lib_filename, 'w') as file:
            file.write(content)
    
    # Clean up
    try:
        import os
        os.unlink(temp_filename)
    except:
        pass
