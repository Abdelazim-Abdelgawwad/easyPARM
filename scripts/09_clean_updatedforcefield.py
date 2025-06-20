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

import periodictable
import re

# Function to read metal numbers from a file
def read_metal_numbers(file_path):
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []

# Function to read new atom type from a file
def read_new_atom_types(file_path):
    try:
        with open(file_path, 'r') as file:
            return {line.strip() for line in file if line.strip()}
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return set()

# Removes lines related to metal atoms from BOND, ANGLE, DIHE, and IMPROPER sections.
# Does not remove any lines from the NONBON section.
def remove_metal_lines(input_file, output_file, metal_atom_types):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_section = None
            lines_kept = 0
            lines_removed = 0

            # Create pattern strings for exact matching of each metal atom type
            # This handles cases like "n1-b" where we need to match 'b' as a complete atom type
            metal_patterns = []
            for metal_type in metal_atom_types:
                # Match the metal type at word boundaries or with hyphens/dots
                pattern = rf'(^|\s|\-|\.)({re.escape(metal_type)})($|\s|\-|\.)'
                metal_patterns.append(pattern)
            
            for line in infile:
                stripped_line = line.strip()

                # Handle section headers
                if stripped_line in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                    current_section = stripped_line
                    outfile.write(line)
                    lines_kept += 1
                    continue

                # Handle empty lines
                if not stripped_line:
                    outfile.write(line)
                    lines_kept += 1
                    continue

                # Check if we should remove this line
                should_remove = False
                
                if current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                    # Use our crafted regex patterns to detect metal atom types
                    for pattern in metal_patterns:
                        if re.search(pattern, stripped_line):
                            should_remove = True
                            matched_metal = re.search(pattern, stripped_line).group(2)
                            lines_removed += 1
                            break
                
                if not should_remove:
                    outfile.write(line)
                    lines_kept += 1

    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")
        raise

# Copy lines related to metal atoms from BOND, ANGLE, DIHE, and IMPROPER sections.
def copy_metal_lines(forcefield_file, output_file, metal_atom_types):
    try:
        with open(forcefield_file, 'r') as infile:
            forcefield_lines = infile.readlines()

        with open(output_file, 'r') as outfile:
            updated_lines = outfile.readlines()

        section_lines = {"BOND": [], "ANGLE": [], "DIHE": [], "IMPROPER": [], "NONBON": []}
        current_section = None
        lines_copied = 0

        # Create pattern strings for exact matching of each metal atom type
        metal_patterns = []
        for metal_type in metal_atom_types:
            # Match the metal type at word boundaries or with hyphens/dots
            pattern = rf'(^|\s|\-|\.)({re.escape(metal_type)})($|\s|\-|\.)'
            metal_patterns.append(pattern)
            
        for line in forcefield_lines:
            stripped_line = line.strip()

            if stripped_line in section_lines.keys():
                current_section = stripped_line
                continue

            if not stripped_line:
                continue

            if current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                # Check if the line contains any of our metal atom types
                contains_metal = False
                for pattern in metal_patterns:
                    if re.search(pattern, stripped_line):
                        contains_metal = True
                        matched_metal = re.search(pattern, stripped_line).group(2)
                        section_lines[current_section].append(line)
                        lines_copied += 1
                        break

        final_lines = []
        current_section = None
        for line in updated_lines:
            stripped_line = line.strip()
            final_lines.append(line)

            if stripped_line in section_lines.keys():
                current_section = stripped_line
                if section_lines[current_section]:
                    final_lines.extend(section_lines[current_section])
                section_lines[current_section] = []

        with open(output_file, 'w') as outfile:
            outfile.writelines(final_lines)

    except Exception as e:
        print(f"Error inserting metal lines into the correct sections: {str(e)}")
        raise

# Extracts the element symbol from an atom name.
# Focuses on the first character(s) that represent the element.
def get_element_from_atom_name(atom_name):
    # Extract the alphabetic prefix (stopping at the first digit)
    match = re.match(r'([A-Za-z]+)', atom_name)
    if match:
        element_name = match.group(1)
        # Try two-letter elements first
        if len(element_name) >= 2:
            possible_element = element_name[:2].capitalize()
            if possible_element in {elem.symbol for elem in periodictable.elements}:
                return possible_element
        # Try one-letter element
        possible_element = element_name[0].upper()
        if possible_element in {elem.symbol for elem in periodictable.elements}:
            return possible_element
    return None

# Creates a mapping of atom types to their corresponding elements based on atom names in mol2 file.
def create_atom_type_to_element_mapping(atom_info):
    atom_type_to_element = {}
    atom_type_to_name = {}
    
    # First, create a mapping of atom types to all their corresponding atom names
    for atom_data in atom_info.values():
        atom_type = atom_data['type']
        atom_name = atom_data['name']
        if atom_type not in atom_type_to_name:
            atom_type_to_name[atom_type] = set()
        atom_type_to_name[atom_type].add(atom_name)
    
    # Then, determine the element for each atom type
    for atom_type, atom_names in atom_type_to_name.items():
        # Get all unique elements from the atom names
        elements = {get_element_from_atom_name(name) for name in atom_names if get_element_from_atom_name(name)}
        # If we get exactly one element type, use that
        if len(elements) == 1:
            element = elements.pop()
            if element:
                atom_type_to_element[atom_type] = element
    
    return atom_type_to_element

# Updates the mass section using atom names from mol2 file to determine correct masses.
# Removes lines for atom types that don't exist in the mol2 file.
def update_mass_section(input_file, output_file, metal_atom_types, atom_info):
    try:
        # Create mapping of atom types to elements based on mol2 atom names
        atom_type_to_element = create_atom_type_to_element_mapping(atom_info)
        
        # Get all unique atom types from mol2 file
        mol2_atom_types = {data['type'] for data in atom_info.values()}

        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            in_mass_section = False
            lines_kept = 0
            lines_removed = 0
            masses_updated = 0
            
            for line in infile:
                stripped_line = line.strip()

                if stripped_line == "MASS":
                    in_mass_section = True
                    outfile.write(line)
                    lines_kept += 1
                    continue

                if in_mass_section:
                    if stripped_line in ["BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                        in_mass_section = False
                        outfile.write(line)
                        lines_kept += 1
                        continue

                    # Skip empty lines
                    if not stripped_line:
                        outfile.write(line)
                        lines_kept += 1
                        continue

                    parts = line.split()
                    if len(parts) >= 3:
                        atom_type = parts[0]
                        
                        # Skip this line if atom type doesn't exist in mol2 file
                        if atom_type not in mol2_atom_types:
                            lines_removed += 1
                            continue
                            
                        try:
                            current_mass = float(parts[1])
                        except ValueError:
                            outfile.write(line)
                            lines_kept += 1
                            continue

                        # Get element symbol from our mapping
                        element_symbol = atom_type_to_element.get(atom_type)
                        
                        if element_symbol:
                            element = getattr(periodictable, element_symbol, None)
                            if element:
                                correct_mass = element.mass
                                if abs(current_mass - correct_mass) > 0.1:
                                    # Preserve formatting while updating mass
                                    mass_start_pos = line.find(str(current_mass))
                                    if mass_start_pos != -1:
                                        mass_end_pos = mass_start_pos + len(str(current_mass))
                                        new_line = (
                                            f"{line[:mass_start_pos]}"
                                            f"{correct_mass:8.3f}"
                                            f"{line[mass_end_pos:]}"
                                        )
                                        outfile.write(new_line)
                                        lines_kept += 1
                                        masses_updated += 1
                                        continue
                        
                    outfile.write(line)
                    lines_kept += 1
                else:
                    outfile.write(line)
                    lines_kept += 1


    except Exception as e:
        print(f"Error updating MASS section: {str(e)}")
        raise

# Extracts atom information from mol2 file.
def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_info = {}
        is_atom_section = False
        atom_count = 0
        
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                is_atom_section = False
                continue
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    full_atom_name = parts[1]
                    atom_type = parts[5]
                    atom_info[atom_id] = {'name': full_atom_name, 'type': atom_type}
                    atom_count += 1
        
        return atom_info
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

# Main function that runs the overall functions
def main():
    try:
        atom_info = extract_atom_types_from_mol2('NEW_COMPLEX.mol2')
        if not atom_info:
            print("Error: No atom information extracted from mol2 file.")
            return

        metal_ids = read_metal_numbers('metal_number.dat')
        if not metal_ids:
            print("Error: No metal atom IDs found.")
            return
        # Get metal atom types
        metal_atom_types = {atom_info[metal_id]['type'] for metal_id in metal_ids if metal_id in atom_info}
        if not metal_atom_types:
            print("Error: Could not determine metal atom types.")
            return
        remove_metal_lines('updated_COMPLEX_modified.frcmod', 'temp_COMPLEX_modified.frcmod', metal_atom_types)
        
        update_mass_section('temp_COMPLEX_modified.frcmod', 'updated_COMPLEX_modified2.frcmod', metal_atom_types, atom_info)
        
        copy_metal_lines('forcefield2.dat', 'updated_COMPLEX_modified2.frcmod', metal_atom_types)
        
        
    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
