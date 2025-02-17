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
#                              |  $$$$$$/              Ver. 3.00 - 12 January 2024                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de CiÃ¨ncia Molecular (ICMol), Universitat de ValÃ¨ncia, P.O. Box 22085, ValÃ¨ncia 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de ValÃ¨ncia. E-mail: abdelazim.abdelgawwad@uv.es         #
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
    "Cl", "Br", "Se", "Ne", "He", "Li", "Mg", "Al",
    "Xe", "Cs", "Ba", "La", "Pr", "Pm", "Sm", "Eu",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Pa",
    "Pu", "Am", "Cf", "Es", "Fm", "Md", "No"
}

# Metal elements (that should remain unchanged)
metal_two_letter_elements = {
    "Ru", "Pd", "Ag", "Pt", "Rh", "Ir", 'Cr', 'Co', 'Re', 'Ir', 'Sn', 'Gd', 'In', 'Sc', 'Ar', 'Fe', 'Zn', 'Si', 'Ni', "Sb", "Ti", "Mn", "Ni", "Cu", "Ga", "Ge" , "As", "Rb", "Sr", "Te", "Au", "Pb", "Hg", "Bi", "Po", "Rn", "Fr", "Ra", "Ac", "Th", "Ta"
}

# Function to generate counters for different metal types
def get_counter_for_metal(metal_index, atom_type):
    if metal_index == 0:
        # Numbers 0-9 followed by lowercase letters a to k
        return iter(list(range(1, 10)) + list('zyxwvu'))
    elif metal_index == 1:
        if atom_type.isupper():
            return iter("abcdefghijk")  # Lowercase letters: a, b, c, ...
        else:
            return iter("ABCDEFGHIJK")  # Uppercase letters: A, B, C, ...
    elif metal_index == 2:
        if atom_type.isupper():
            # Forward alphabetical lowercase: l through v
            return iter("lmnopqrstuv")
        else:
            # Forward alphabetical uppercase: L through V
            return iter("LMNOPQRSTUV")
    elif metal_index == 4:
        # Numbers 9-0 followed by lowercase letters n to z
            return iter(list(range(9, -1, -1)) + list('nopqrstuvwxyz'))
    elif metal_index == 3:
        if atom_type.isupper():
            # Sequence starting from middle: n through z, then a through m
            return iter("nopqrstuvwxyzabcdefghijklm")
        else:
            # Uppercase sequence starting from middle: N through Z, then A through M
            return iter("NOPQRSTUVWXYZABCDEFGHIJKLM")
    elif metal_index == 5:
        if atom_type.isupper():
            # Forward alphabetical lowercase: w through z, then a through m
            return iter("wxyzabcdefghijklm")
        else:
            # Forward alphabetical uppercase: W through Z, then A through M
            return iter("WXYZABCDEFGHIJKLM")
    else:
        raise ValueError(f"No counter available for metal index {metal_index}")
    
# Function to update the mol2 file with new atom types
def update_mol2_file(atom_positions, mol2_file, new_file, metal_positions, bonds):
    with open(mol2_file, 'r') as file:
        lines = file.readlines()
    # Find the start of atom and bond sections
    atom_start = lines.index("@<TRIPOS>ATOM\n")
    bond_start = lines.index("@<TRIPOS>BOND\n")
    atom_lines = lines[atom_start + 1:bond_start]
    metal_counters = {}
    two_letter_counters = {}
    global_counter = 1
    updated_lines = []
    new_atom_types = []
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
                    if normalized_atom_type in normal_two_letter_elements:
                        # Use or initialize a counter for this specific non-metal two-letter element
                        if normalized_atom_type not in two_letter_counters:
                            two_letter_counters[normalized_atom_type] = 1
                        new_atom_type = f"{normalized_atom_type}{two_letter_counters[normalized_atom_type]}"
                        two_letter_counters[normalized_atom_type] += 1
                    elif normalized_atom_type in metal_two_letter_elements:
                        # Keep the metal element unchanged
                        new_atom_type = normalized_atom_type
                    else:
                        # Initialize counter for this metal and atom type
                        if metal_index not in metal_counters:
                            metal_counters[metal_index] = get_counter_for_metal(metal_index, atom_type)
                        # Generate new atom type based on original type
                        if atom_type.isupper():
                            new_atom_type = f"{atom_type[0]}{next(metal_counters[metal_index])}".lower()
                        else:
                            new_atom_type = f"{atom_type[0]}{next(metal_counters[metal_index])}".upper()
                    atom_type = new_atom_type
                    new_atom_types.append(new_atom_type)
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
    
    # Write new atom types to new_atomtype.dat
    with open('new_atomtype.dat', 'w') as file:
        for new_atom_type in new_atom_types:
            file.write(new_atom_type + '\n')
    
    # Write hybridization information to Hybridization_info.dat
    with open('Hybridization_Info.dat', 'w') as file:
        for atom_name, new_atom_type in zip(atom_names, new_atom_types):
            file.write(f"{new_atom_type}  {atom_name}  sp3\n")
# Main function to excute the process
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


