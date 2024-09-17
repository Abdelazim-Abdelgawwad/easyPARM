import os
import re

# Extract the info from the MOL2 file 
def extract_atom_info_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_info = {}
        is_atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    # Remove any numbers from the atom name
                    atom_name_base = re.sub(r'\d+', '', atom_name)
                    atom_info[atom_id] = atom_name_base
        return atom_info
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

# Read the metal number to detect the bonded atom
def read_metal_numbers(file_path):
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []

# Calculate the number of bonds linked 
def count_bonds(metal_numbers, distance_file_path):
    bond_counts = {metal: 0 for metal in metal_numbers}
    try:
        with open(distance_file_path, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) >= 2:
                    atom1 = int(parts[0])
                    atom2 = int(parts[1])
                    if atom1 in bond_counts:
                        bond_counts[atom1] += 1
                    if atom2 in bond_counts:
                        bond_counts[atom2] += 1
        return bond_counts
    except Exception as e:
        print(f"Error reading {distance_file_path}: {str(e)}")
        return {}

# Read UFF data from the existing file
def read_uff_data():
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        uff_data_file_path = os.path.join(script_dir, 'uff_data.txt')
        
        with open(uff_data_file_path, 'r') as file:
            uff_data = {}
            for line in file:
                parts = line.split()
                if len(parts) >= 7:
                    atom_type = parts[0]
                    col4_half = float(parts[3]) / 2.0
                    col5 = float(parts[4])
                    uff_data[atom_type] = (col4_half, col5)
        return uff_data
    except Exception as e:
        print(f"Error reading {uff_data_file_path}: {str(e)}")
        return {}

# Update the NONBON section with the correct info from UFF file
def update_frcmod_with_uff_data(frcmod_file_path, uff_data, atom_info, bond_counts):
    try:
        with open(frcmod_file_path, 'r') as file:
            lines = file.readlines()

        nonbon_start_index = None
        for i, line in enumerate(lines):
            if line.strip().startswith("NONBON"):
                nonbon_start_index = i + 1
                break

        if nonbon_start_index is None:
            print("No NONBON section found in the frcmod file.")
            return

        for i in range(nonbon_start_index, len(lines)):
            parts = lines[i].split()
            if len(parts) >= 3:
                old_atom_type = parts[0].strip()
                matching_atom_nums = [num for num, aname in atom_info.items() if aname == old_atom_type]
                if not matching_atom_nums:
                    continue
                
                atom_num = matching_atom_nums[0]
                if atom_num in bond_counts:
                    bond_count = bond_counts[atom_num]
                    base_atom_type = old_atom_type

                    # Find the correct new atom type by looking for a match in uff_data
                    matching_uff_keys = [key for key in uff_data.keys() if key.startswith(base_atom_type)]
                    new_atom_type = None

                    # Match the bond count + suffix
                    for key in matching_uff_keys:
                        if key[len(base_atom_type):].startswith(str(bond_count)):
                            new_atom_type = key
                            break

                    if new_atom_type and new_atom_type in uff_data:
                        col2, col3 = uff_data[new_atom_type]
                    else:
                        # Fallback to the base atom type only
                        matched_key = None
                        for key in uff_data.keys():
                            if key.startswith(base_atom_type) and key != new_atom_type:
                                matched_key = key
                                break
                        
                        if matched_key:
                            col2, col3 = uff_data[matched_key]
                        else:
                            continue

                    # Update the frcmod line with the new values
                    lines[i] = f"{old_atom_type:<10}{col2:>8.4f} {col3:>8.4f}\n"

        with open(frcmod_file_path, 'w') as file:
            file.writelines(lines)

    except Exception as e:
        print(f"Error processing {frcmod_file_path}: {str(e)}")

# Inputs
mol2_file_path = 'NEW_COMPLEX.mol2'
metal_numbers_file_path = 'metal_number.dat'
distance_file_path = 'distance.dat'
frcmod_file_path = 'updated_updated_COMPLEX_modified2.frcmod'

# Run the extraction of info and total bond number
atom_info = extract_atom_info_from_mol2(mol2_file_path)
metal_numbers = read_metal_numbers(metal_numbers_file_path)
bond_counts = count_bonds(metal_numbers, distance_file_path)

# Update the file
uff_data = read_uff_data()
update_frcmod_with_uff_data(frcmod_file_path, uff_data, atom_info, bond_counts)
