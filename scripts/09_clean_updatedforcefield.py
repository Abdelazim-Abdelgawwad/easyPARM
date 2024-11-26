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
#                              |  $$$$$$/              Ver. 2.00 - 1 November 2024                                #
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

#Removes lines related to metal atoms from BOND, ANGLE, DIHE, and IMPROPER sections.
#Does not remove any lines from the NONBON section.
def remove_metal_lines(input_file, output_file, metal_atom_types):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_section = None

            for line in infile:
                stripped_line = line.strip()

                if stripped_line in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                    current_section = stripped_line
                    outfile.write(line)
                    continue

                if not stripped_line:
                    outfile.write(line)
                    continue

                atom_types_in_line = re.findall(r'\b[A-Z][a-zA-Z0-9]*\b', stripped_line)

                if current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                    if any(atom_type in metal_atom_types for atom_type in atom_types_in_line):
                        continue

                outfile.write(line)

    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")

#Copy lines related to metal atoms from BOND, ANGLE, DIHE, and IMPROPER sections.
def copy_metal_lines(forcefield_file, output_file, metal_atom_types):
    try:
        with open(forcefield_file, 'r') as infile:
            forcefield_lines = infile.readlines()

        with open(output_file, 'r') as outfile:
            updated_lines = outfile.readlines()

        section_lines = {"BOND": [], "ANGLE": [], "DIHE": [], "IMPROPER": [], "NONBON": []}
        current_section = None

        for line in forcefield_lines:
            stripped_line = line.strip()

            if stripped_line in section_lines.keys():
                current_section = stripped_line
                continue

            if current_section in ["BOND", "ANGLE", "DIHE", "IMPROPER"]:
                if any(atom_type in stripped_line for atom_type in metal_atom_types):
                    section_lines[current_section].append(line)

        final_lines = []
        current_section = None
        for line in updated_lines:
            stripped_line = line.strip()
            final_lines.append(line)

            if stripped_line in section_lines.keys():
                current_section = stripped_line
                final_lines.extend(section_lines[current_section])
                section_lines[current_section] = []

        with open(output_file, 'w') as outfile:
            outfile.writelines(final_lines)

    except Exception as e:
        print(f"Error inserting metal lines into the correct sections: {str(e)}")

#Extracts the element symbol from an atom name.
#Focuses on the first character(s) that represent the element.
def get_element_from_atom_name(atom_name):
    # Extract the alphabetic prefix (stopping at the first digit)
    match = re.match(r'([A-Za-z]+)', atom_name)
    if match:
        element_name = match.group(1)
        if len(element_name) >= 2 and element_name[:2] in {elem.symbol for elem in periodictable.elements}:
            return element_name[:2]
        return element_name[0]
    return None

#Creates a mapping of atom types to their corresponding elements based on atom names in mol2 file.
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
        elements = {get_element_from_atom_name(name) for name in atom_names}
        # If we get exactly one element type, use that
        if len(elements) == 1:
            element = elements.pop()
            if element:
                atom_type_to_element[atom_type] = element
    
    return atom_type_to_element

#Updates the mass section using atom names from mol2 file to determine correct masses.
#Removes lines for atom types that don't exist in the mol2 file.
def update_mass_section(input_file, output_file, metal_atom_types, atom_info):
    try:
        # Create mapping of atom types to elements based on mol2 atom names
        atom_type_to_element = create_atom_type_to_element_mapping(atom_info)
        
        # Get all unique atom types from mol2 file
        mol2_atom_types = {data['type'] for data in atom_info.values()}

        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            in_mass_section = False
            for line in infile:
                stripped_line = line.strip()

                if stripped_line == "MASS":
                    in_mass_section = True
                    outfile.write(line)
                    continue

                if in_mass_section:
                    if stripped_line in ["BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                        in_mass_section = False
                        outfile.write(line)
                        continue

                    # Skip empty lines
                    if not stripped_line:
                        outfile.write(line)
                        continue

                    parts = line.split()
                    if len(parts) >= 3:
                        atom_type = parts[0]
                        
                        # Skip this line if atom type doesn't exist in mol2 file
                        if atom_type not in mol2_atom_types:
                            continue
                            
                        try:
                            current_mass = float(parts[1])
                        except ValueError:
                            outfile.write(line)
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
                                        continue
                        
                    outfile.write(line)
                else:
                    outfile.write(line)

    except Exception as e:
        print(f"Error updating MASS section: {str(e)}")
        raise

#Extracts atom information from mol2 file.
def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_info = {}
        is_atom_section = False
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
        return atom_info
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

#Main function that run the overall functions
def main():
    atom_info = extract_atom_types_from_mol2('NEW_COMPLEX.mol2')

    metal_ids = read_metal_numbers('metal_number.dat')
    metal_atom_types = {atom_info[metal_id]['type'] for metal_id in metal_ids if metal_id in atom_info}

    remove_metal_lines('updated_COMPLEX_modified.frcmod', 'temp_COMPLEX_modified.frcmod', metal_atom_types)
    update_mass_section('temp_COMPLEX_modified.frcmod', 'updated_COMPLEX_modified2.frcmod', metal_atom_types, atom_info)
    copy_metal_lines('forcefield2.dat', 'updated_COMPLEX_modified2.frcmod', metal_atom_types)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()
