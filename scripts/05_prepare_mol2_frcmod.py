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

# Function to read metal numbers from a file
def read_metal_numbers(file_path):
    with open(file_path, 'r') as file:
        # Read each line, strip whitespace, and convert to integer
        metal_numbers = [int(line.strip()) for line in file]
    return metal_numbers

# Function to read bond distances and create a dictionary of bonds
def read_distances(file_path):
    bonds = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            pos1, pos2 = int(parts[0]), int(parts[1])
            # Create bidirectional bonds
            if pos1 not in bonds:
                bonds[pos1] = []
            if pos2 not in bonds:
                bonds[pos2] = []
            bonds[pos1].append(pos2)
            bonds[pos2].append(pos1)
    return bonds

# Function to extract atoms bonded to metal positions
def extract_atoms(metal_positions, bonds):
    bonded_atoms = set()
    for metal in metal_positions:
        # Add all atoms bonded to each metal to the set
        bonded_atoms.update(bonds.get(metal, []))
    return bonded_atoms

# List of two-letter elements from the periodic table, using uppercase for consistency
normal_two_letter_elements = {
    "Xe"#"Cl", "Br", "Se", "Ne", "He", "Li", "Mg", "Al",
    #"Xe", "Cs", "Ba", "La", "Pr", "Pm", "Sm", "Eu",
    #"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    #"Pa", "Pu", "Am", "Cf", "Es", "Fm", "Md", "No"
}

# Metal elements (that should remain unchanged)
metal_two_letter_elements = {
    "Ru", "Pd", "Ag", "Pt", "Rh", "Zr", "Ir", 'Cr', 'Co', 'Re', 'Ir', 'Sn', 'Gd', 'In', 'Sc', 'Ar', 'Fe', 'Zn', 'Si', 'Ni', "Sb", "Ti", "Mn", "Cu", "Ga", "Ge" , "As", "Rb", "Sr", "Te", "Au", "Pb", "Hg", "Bi", "Po", "Rn", "Fr", "Ra", "Ac", "Th", "Ta"
}

# Keep track of all new_atom_types created
used_atom_types = set()

# Create a master list of all available counters for fallback
all_available_counters = list('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')

# Function to read existing atom types from the mol2 file
def read_existing_atom_types(mol2_file):
    existing_types = set()
    with open(mol2_file, 'r') as file:
        lines = file.readlines()
    
    # Find the atom section
    try:
        atom_start = lines.index("@<TRIPOS>ATOM\n")
        bond_start = lines.index("@<TRIPOS>BOND\n")
        atom_lines = lines[atom_start + 1:bond_start]
        
        # Extract atom types from each line
        for line in atom_lines:
            parts = line.split()
            if len(parts) >= 6:  # Ensure there are enough columns
                atom_type = parts[5]
                existing_types.add(atom_type)
    except ValueError:
        # Handle case where the file format is unexpected
        print("Warning: Could not parse mol2 file sections properly")
    
    return existing_types

#Generate a unique sequence of 12 single-character counters for each metal index.
#Each metal index gets its own set of characters (with some reuse across different indices).
def get_counter_for_metal(metal_index, atom_type):
    # Dividing standard alphanumeric characters among 12 metal indices
    if metal_index == 1:
        # Digits 0-9 + a,b (12 characters)
        return iter(list('1234567891ab'))
    elif metal_index == 2:
        # Uppercase A through L (12 characters)
        return iter(list('ABCDEFGHIJKL'))
    elif metal_index == 3:
        # Lowercase c through n (12 characters)
        return iter(list('abcdefghijklmn'))
    elif metal_index == 4:
        # Lowercase o through z (12 characters)
        return iter(list('opqrstuvwxyz'))
    elif metal_index == 5:
        # Uppercase M through X (12 characters) 
        return iter(list('MNOPQRSTUVWX'))
    elif metal_index == 6:
        # Remaining uppercase Y,Z + additional digits and letters to reach 12
        return iter(list('YZabcdefghij'))
    elif metal_index == 7:
        # Another set reusing alphanumeric chars
        return iter(list('klmnopqrstuv'))
    elif metal_index == 8:
        # Another set reusing alphanumeric chars
        return iter(list('wxyzABCDEFGH'))
    elif metal_index == 9:
        # Another set reusing alphanumeric chars
        return iter(list('IJKLMNOPQRST'))
    elif metal_index == 10:
        # Another set reusing alphanumeric chars
        return iter(list('UVWXYZ012345'))
    elif metal_index == 11:
        # Another set reusing alphanumeric chars
        return iter(list('6789abcdefgh'))
    elif metal_index == 12:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 13:
        # Remaining uppercase Y,Z + additional digits and letters to reach 12
        return iter(list('YZabcdefghij'))
    elif metal_index == 14:
        # Another set reusing alphanumeric chars
        return iter(list('klmnopqrstuv'))
    elif metal_index == 15:
        # Another set reusing alphanumeric chars
        return iter(list('wxyzABCDEFGH'))
    elif metal_index == 16:
        # Another set reusing alphanumeric chars
        return iter(list('IJKLMNOPQRST'))
    elif metal_index == 17:
        # Another set reusing alphanumeric chars
        return iter(list('UVWXYZ012345'))
    elif metal_index == 18:
        # Another set reusing alphanumeric chars
        return iter(list('6789abcdefgh'))
    elif metal_index == 19:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 20:
        # Another set reusing alphanumeric chars
        return iter(list('6789abcdefgh'))
    elif metal_index == 22:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 23:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 24:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 25:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    elif metal_index == 26:
        # Another set reusing alphanumeric chars
        return iter(list('ijklmnopqrst'))
    else:
        # Fallback for any other index
        raise ValueError(f"No counter available for metal index {metal_index}")

#Generate a unique new_atom_type using counters, ensuring no duplicates.
#If the first choice is already used, find another unused letter from the available sets.
def get_unique_new_atom_type(atom_type, metal_index):
    # Get the primary counters for this metal index
    counter_iter = get_counter_for_metal(metal_index, atom_type)
    
    # Try each counter in sequence until we find one that creates a unique new_atom_type
    for counter in counter_iter:
        new_atom_type = f"{atom_type[0]}{counter}"
        
        if new_atom_type not in used_atom_types:
            # Found a unique new_atom_type, add it to the set and return it
            used_atom_types.add(new_atom_type)
            return new_atom_type
    
    # If we've used all counters from the preferred set for this metal,
    # let's find any unused counter from our master list
    for counter in all_available_counters:
        new_atom_type = f"{atom_type[0]}{counter}"
        
        if new_atom_type not in used_atom_types:
            used_atom_types.add(new_atom_type)
            return new_atom_type
    
    # If even the master list is exhausted, only then resort to numeric suffix
    i = 1
    while True:
        new_atom_type = f"{atom_type[0]}{i}"
        if new_atom_type not in used_atom_types:
            used_atom_types.add(new_atom_type)
            return new_atom_type
        i += 1

# Function to update the mol2 file with new atom types
def update_mol2_file(atom_positions, mol2_file, new_file, metal_positions, bonds):
    # First, read all existing atom types from the mol2 file
    existing_types = read_existing_atom_types(mol2_file)
    
    # Add existing types to our used_atom_types set
    global used_atom_types
    used_atom_types.update(existing_types)
    
    with open(mol2_file, 'r') as file:
        lines = file.readlines()
    # Find the start of atom and bond sections
    atom_start = lines.index("@<TRIPOS>ATOM\n")
    bond_start = lines.index("@<TRIPOS>BOND\n")
    atom_lines = lines[atom_start + 1:bond_start]
    two_letter_counters = {}
    updated_lines = []
    new_atom_types = []
    original_atom_types = []  # New list to store original atom types
    atom_names = []  # New list to store atom names
    
    for line in atom_lines:
        parts = line.split()
        # Extract columns from the line
        atom_id = int(parts[0])
        atom_name, x, y, z, atom_type, int_number, name, charge = parts[1:9]
        x, y, z, charge = map(float, (x, y, z, charge))
        int_number = int(int_number)
        # Normalize the atom type (for comparison purposes)
        normalized_atom_type = atom_type.capitalize()
        # Check if the atom needs to be modified
        if atom_id in atom_positions:
            for metal_index, metal in enumerate(metal_positions):
                if atom_id in bonds.get(metal, []):
                    # Store the original atom type before any modification
                    original_atom_type = atom_type
                    
                    if normalized_atom_type in normal_two_letter_elements:
                        # Use or initialize a counter for this specific non-metal two-letter element
                        if normalized_atom_type not in two_letter_counters:
                            two_letter_counters[normalized_atom_type] = 1
                        new_atom_type = f"{normalized_atom_type}{two_letter_counters[normalized_atom_type]}"
                        
                        # Check if this new_atom_type already exists in our comprehensive list
                        while new_atom_type in used_atom_types:
                            two_letter_counters[normalized_atom_type] += 1
                            new_atom_type = f"{normalized_atom_type}{two_letter_counters[normalized_atom_type]}"
                            
                        two_letter_counters[normalized_atom_type] += 1
                    elif normalized_atom_type in metal_two_letter_elements:
                        # Keep the metal element unchanged
                        new_atom_type = normalized_atom_type
                    else:
                        # Use the new function to generate a unique atom type
                        # metal_index+1 to match the new function's expected range (1-12)
                        new_atom_type = get_unique_new_atom_type(atom_type, metal_index+1)
                    
                    atom_type = new_atom_type
                    new_atom_types.append(new_atom_type)
                    original_atom_types.append(original_atom_type)  # Store the original atom type
                    # Store just the letter part of the atom name (remove numbers)
                    base_atom_name = ''.join(c for c in atom_name if not c.isdigit())
                    atom_names.append(base_atom_name)  # Store the corresponding atom name without numbers
                    break
        # Reformat the line preserving the original structure
        updated_line = f"{atom_id:7d} {atom_name:<4s} {x:10.4f} {y:10.4f} {z:10.4f} {atom_type:<6s} {int_number:3d} {name:<4s} {charge:10.6f}\n"
        updated_lines.append(updated_line)
    
    # Write the updated content to the new file
    with open(new_file, 'w') as file:
        file.writelines(lines[:atom_start + 1])  # Write everything before atoms
        file.writelines(updated_lines)  # Write modified atoms
        file.writelines(lines[bond_start:])  # Write everything after atoms
    
    # Write new atom types and original atom types to new_atomtype.dat
    with open('new_atomtype.dat', 'w') as file:
        for new_type, orig_type in zip(new_atom_types, original_atom_types):
            file.write(f"{new_type} \n")
    
    # Write new atom types and original atom types to new_atomtype.dat
    with open('metalloprotein_atomtype.dat', 'w') as file:
        for new_type, orig_type in zip(new_atom_types, original_atom_types):
            file.write(f"{new_type} {orig_type}\n")
    
    # Write hybridization information to Hybridization_info.dat
    with open('Hybridization_Info.dat', 'w') as file:
        file.write("addAtomTypes {\n")
        for atom_name, new_atom_type in zip(atom_names, new_atom_types):
            file.write(f'  {{ "{new_atom_type}" "{atom_name}" "sp3" }}\n')
        file.write("}\n")

# Main function to execute the process
def main():
    metal_numbers_file = 'metal_number.dat'
    distances_file = 'distance.dat'
    mol2_file = 'COMPLEX_modified.mol2'
    new_file = 'NEW_COMPLEX.mol2'

    # Read input files
    metal_positions = read_metal_numbers(metal_numbers_file)
    bonds = read_distances(distances_file)

    # Extract atoms bonded to metals
    atom_positions = extract_atoms(metal_positions, bonds)
    
    # Update the mol2 file and save it as a new file
    update_mol2_file(atom_positions, mol2_file, new_file, metal_positions, bonds)

if __name__ == "__main__":
    main()
