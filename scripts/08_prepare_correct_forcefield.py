import itertools
import re

# Function to extract atom types from a MOL2 file
def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_types = {}
        is_atom_section = False
        for line in lines:
            # Check for the start of the atom section
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            # Check for the end of the atom section
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_type = parts[5]
                    atom_types[atom_id] = atom_type
        return atom_types
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

# Function to read metal numbers from a file
def read_metal_numbers(file_path):
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []

# Function to read distance information from a file
def read_distance_file(file_path):
    bonds = {}
    try:
        with open(file_path, 'r') as file:
            file_contents = file.read()
            
            for line in file_contents.split('\n'):
                parts = line.split()
                if len(parts) == 3:
                    atom1, atom2, distance = map(float, parts)
                    bonds[int(atom1), int(atom2)] = distance
        return bonds
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

# Function to update the FRCMOD file with new atom types
def update_frcmod_file(input_file, output_file, old_to_new_types, metal_specific_types):
    def match_atom_type(atom_type, pattern):
        # Convert pattern to regex, replacing '*' with '.*'
        regex_pattern = '^' + re.escape(pattern).replace(r'\*', '.*') + '$'
        return re.match(regex_pattern, atom_type, re.IGNORECASE) is not None

    def replace_atom_type(line, old_type, new_type, current_section):
        # Find the position of the atom type (case-insensitive)
        match = re.search(re.escape(old_type), line, re.IGNORECASE)
        if match:
            index = match.start()
            # Preserve the original format by replacing only the atom type, maintaining all spaces
            return line[:index] + new_type + line[index + len(old_type):]
        return line

    def generate_replacements(atoms, old_to_new_types, current_section):
        replacements = []
        replaceable_indices = [i for i, atom in enumerate(atoms) if any(match_atom_type(atom, old_type) for old_type in old_to_new_types)]
        
        # Generate single replacements
        for i in replaceable_indices:
            for old_type, new_types in old_to_new_types.items():
                if match_atom_type(atoms[i], old_type):
                    for new_type in new_types:
                        new_atoms = atoms.copy()
                        new_atoms[i] = new_type
                        replacements.append(new_atoms)
        
        # Generate combinations of replacements only for DIHE section
        if current_section == "DIHE":
            for r in range(2, len(replaceable_indices) + 1):
                for indices in itertools.combinations(replaceable_indices, r):
                    combinations = []
                    for idx in indices:
                        atom_combinations = []
                        for old_type, new_types in old_to_new_types.items():
                            if match_atom_type(atoms[idx], old_type):
                                atom_combinations.extend(new_types)
                        combinations.append(atom_combinations)
                    
                    for combination in itertools.product(*combinations):
                        new_atoms = atoms.copy()
                        for idx, new_type in zip(indices, combination):
                            new_atoms[idx] = new_type
                        replacements.append(new_atoms)
        
        return replacements

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_section = None
            for line in infile:
                original_line = line.rstrip()

                if original_line in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                    current_section = original_line
                    outfile.write(original_line + '\n')
                    continue

                if not original_line or original_line.startswith('Remark') or original_line.startswith('#'):
                    outfile.write(original_line + '\n')
                    continue

                parts = original_line.split()
                if not parts:
                    outfile.write(original_line + '\n')
                    continue

                if current_section in ["MASS", "NONBON"]:
                    atom_type = parts[0]
                    outfile.write(original_line + '\n')  # Keep the original line
                    for old_type, new_types in old_to_new_types.items():
                        if match_atom_type(atom_type, old_type):
                            for new_type in new_types:
                                new_line = replace_atom_type(original_line, atom_type, new_type, current_section)
                                if new_line != original_line:
                                    outfile.write(new_line + '\n')

                elif current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                    atoms = parts[0].split('-')
                    outfile.write(original_line + '\n')  # Keep the original line
                    replacements = generate_replacements(atoms, old_to_new_types, current_section)
                    for new_atoms in replacements:
                        new_atoms_string = '-'.join(new_atoms)
                        new_line = replace_atom_type(original_line, parts[0], new_atoms_string, current_section)
                        if new_line != original_line:
                            outfile.write(new_line + '\n')
                else:
                    outfile.write(original_line + '\n')

    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")

def main():
    # Read atom types from both MOL2 files
    old_atom_types = extract_atom_types_from_mol2('COMPLEX.mol2')
    new_atom_types = extract_atom_types_from_mol2('NEW_COMPLEX.mol2')
    
    # Read metal numbers
    metal_ids = read_metal_numbers('metal_number.dat')

    # Read distance file
    bonds = read_distance_file('distance.dat')

    # Find atoms linked to each metal
    metal_bonded_atoms = {metal_id: set() for metal_id in metal_ids}
    for (atom1, atom2), distance in bonds.items():
        if atom1 in metal_ids:
            metal_bonded_atoms[atom1].add(atom2)
        if atom2 in metal_ids:
            metal_bonded_atoms[atom2].add(atom1)
    
    # Create mapping of old to new atom types for atoms bonded to metals
    old_to_new_types = {}
    metal_specific_types = {metal_id: set() for metal_id in metal_ids}
    
    for metal_id, bonded_atoms in metal_bonded_atoms.items():
        for atom_id in bonded_atoms:
            if atom_id in old_atom_types and atom_id in new_atom_types:
                old_type = old_atom_types[atom_id]
                new_type = new_atom_types[atom_id]
                if old_type != new_type:
                    if old_type not in old_to_new_types:
                        old_to_new_types[old_type] = set()
                    old_to_new_types[old_type].add(new_type)
                    metal_specific_types[metal_id].add(new_type)

    # Update FRCMOD file
    update_frcmod_file('COMPLEX_modified.frcmod', 'updated_COMPLEX_modified.frcmod', old_to_new_types, metal_specific_types)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()

