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


# Function to read metal numbers from a file
def read_metal_numbers(file_path):
    with open(file_path, 'r') as file:
        metal_numbers = [int(line.strip()) for line in file]
    return metal_numbers

# Function to read bond distances and create a dictionary of bonds
def read_distances(file_path):
    bonds = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            pos1, pos2 = int(parts[0]), int(parts[1])
            if pos1 not in bonds:
                bonds[pos1] = []
            if pos2 not in bonds:
                bonds[pos2] = []
            bonds[pos1].append(pos2)
            bonds[pos2].append(pos1)
    return bonds

#Find groups of metals that are connected to each other
def find_connected_metals(metal_positions, bonds):
    metal_groups = []
    visited = set()

    def dfs(metal, current_group):
        visited.add(metal)
        current_group.append(metal)
        for neighbor in bonds.get(metal, []):
            if neighbor in metal_positions and neighbor not in visited:
                dfs(neighbor, current_group)

    for metal in metal_positions:
        if metal not in visited:
            current_group = []
            dfs(metal, current_group)
            metal_groups.append(sorted(current_group))

    return metal_groups

#Create a mapping of metal positions to their new atom types
def create_metal_mapping(metal_groups, mol2_file):
    # Read current atom types for metals
    current_types = {}
    with open(mol2_file, 'r') as file:
        lines = file.readlines()
        atom_start = lines.index("@<TRIPOS>ATOM\n")
        bond_start = lines.index("@<TRIPOS>BOND\n")
        atom_lines = lines[atom_start + 1:bond_start]

        for line in atom_lines:
            parts = line.split()
            atom_id = int(parts[0])
            atom_type = parts[5]
            current_types[atom_id] = atom_type

    # Create mapping for connected metals
    metal_mapping = {}
    letter_counters = {}  # Counter for each single-letter metal type

    for group in metal_groups:
        for metal in group:
            current_type = current_types[metal]
            if current_type in metal_two_letter_elements:
                # Keep two-letter metals unchanged
                metal_mapping[metal] = current_type
            else:
                # For single-letter metals, append lowercase letters
                base_type = current_type[0]  # Get first letter
                if base_type not in letter_counters:
                    letter_counters[base_type] = iter('abcdefghijklmnopqsturvwxyz')
                metal_mapping[metal] = f"{base_type}{next(letter_counters[base_type])}"

    return metal_mapping

# List of two-letter elements from the periodic table, using uppercase for consistency
normal_two_letter_elements = {
    "Cl", "Br", "Se", "Ne", "He", "Li", "Mg", "Al",
    "Xe", "Cs", "Ba", "La", "Pr", "Pm", "Sm", "Eu",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "Pa", "Pu", "Am", "Cf", "Es", "Fm", "Md", "No"
}

# Metal elements (that should remain unchanged)
metal_two_letter_elements = {
    "Ru", "Pd", "Ag", "Pt", "Rh", "Zr", "Ir", 'Cr', 'Co', 'Re', 'Ir', 'Sn', 'Gd', 'In', 'Sc', 'Ar', 'Fe', 'Zn', 'Si', 'Ni', "Sb", "Ti", "Mn", "Cu", "Ga", "Ge" , "As", "Rb", "Sr", "Te", "Au", "Pb", "Hg", "Bi", "Po", "Rn", "Fr", "Ra", "Ac", "Th", "Ta"
}

# Function to update the mol2 file with new atom types
def update_mol2_file(metal_mapping, mol2_file, new_file, metal_positions):
    with open(mol2_file, 'r') as file:
        lines = file.readlines()

    atom_start = lines.index("@<TRIPOS>ATOM\n")
    bond_start = lines.index("@<TRIPOS>BOND\n")
    atom_lines = lines[atom_start + 1:bond_start]

    updated_lines = []
    new_atom_types = []

    for line in atom_lines:
        parts = line.split()
        atom_id = int(parts[0])
        atom_name, x, y, z, atom_type, int_number, name, charge = parts[1:9]
        x, y, z, charge = map(float, (x, y, z, charge))
        int_number = int(int_number)

        # Only modify atom type if it's a metal position
        if atom_id in metal_mapping:
            atom_type = metal_mapping[atom_id]
            new_atom_types.append(atom_type)

        updated_line = f"{atom_id:7d} {atom_name:<4s} {x:10.4f} {y:10.4f} {z:10.4f} {atom_type:<6s} {int_number:3d} {name:<4s} {charge:10.6f}\n"
        updated_lines.append(updated_line)

    # Write the updated content
    with open(new_file, 'w') as file:
        file.writelines(lines[:atom_start + 1])
        file.writelines(updated_lines)
        file.writelines(lines[bond_start:])

    # Write new atom types
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

    # Find connected metal groups
    metal_groups = find_connected_metals(metal_positions, bonds)

    # Create mapping for metal renumbering
    metal_mapping = create_metal_mapping(metal_groups, mol2_file)

    # Update the mol2 file
    update_mol2_file(metal_mapping, mol2_file, new_file, metal_positions)

if __name__ == "__main__":
    main()
