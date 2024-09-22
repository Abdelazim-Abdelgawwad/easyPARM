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
    "Cl", "Br", "Se", "Ne", "He", "Li", "Be", "Mg", "Al",
    "Ti", "Mn", "Ni", "Cu", "Ga", "Ge", "As", "Kr", "Rb", "Sr",
    "Sb", "Te", "Xe", "Cs", "Ba", "La", "Pr", "Pm", "Sm", "Eu",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    "Pu", "Am", "Bk", "Cf", "Es", "Fm", "Md", "No"
}

# Metal elements (that should remain unchanged)
metal_two_letter_elements = {
    "Ru", "Pd", "Ag", "Pt", "Rh", "Ir", 'Cr', 'Co', 'Re', 'Ir', 'Sn', 'Gd', 'In', 'Sc', 'Ar', 'Fe', 'Zn', 'Si', 'Ni'
}

# Function to generate counters for different metal types
def get_counter_for_metal(metal_index, atom_type):
    if metal_index == 0:
        return iter(range(1, 100))  # Numbers: 1, 2, 3, ...
    elif metal_index == 1:
        if atom_type.isupper():
            return iter("abdefghijklmnopqrstuvwxyz")  # Lowercase letters: a, b, c, ...
        else:
            return iter("ABDEFGHIJKLMNOPQRSTUVWXYZ")  # Uppercase letters: A, B, C, ...
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

