import periodictable
import re

#Extract atom types and names from a MOL2 file.
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
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_name = ''.join([c for c in parts[1] if not c.isdigit()])  # Remove numbers
                    atom_type = parts[5]
                    atom_info[atom_id] = {'name': atom_name, 'type': atom_type}
        return atom_info
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

#Reads the metal numbers from a file.
def read_metal_numbers(file_path):
    
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []

#Reads new atom types from a file.
def read_new_atom_types(file_path):
    
    try:
        with open(file_path, 'r') as file:
            return {line.strip() for line in file if line.strip()}
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return set()

# Function to do
#1- Removes lines related to metal atoms from BOND, ANGLE, DIHE, and IMPROPER sections.
#2- Removes new atom type lines only from the ANGLE section.
#3- Does not remove any lines from the NONBON section.
#4- Handles cases with various formats and additional information.
def remove_metal_and_new_atom_type_lines(input_file, output_file, metal_atom_types, new_atom_types):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_section = None

            for line in infile:
                stripped_line = line.strip()

                # Identify section headers
                if stripped_line in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                    current_section = stripped_line
                    outfile.write(line)
                    continue

                # Skip empty lines
                if not stripped_line:
                    outfile.write(line)
                    continue

                if current_section == "ANGLE":
                    # Extract atom types from the ANGLE line
                    angle_types = stripped_line.split()[0].split('-')
                    
                    # Check if any of the atom types in the ANGLE line are metal or new atom types
                    if any(at.lower() in (mt.lower() for mt in metal_atom_types) or 
                           at.lower() in (nt.lower() for nt in new_atom_types) 
                           for at in angle_types):
                        continue

                elif current_section in ["BOND", "DIHE", "IMPROPER"]:
                    # Extract atom types from the line for other sections
                    atom_types_in_line = re.findall(r'\b[A-Z][a-z]?\d*\b', stripped_line)
                    
                    # Remove lines containing metal atom types
                    if any(at.lower() in (mt.lower() for mt in metal_atom_types) for at in atom_types_in_line):
                        continue

                # For NONBON and MASS sections, don't remove any lines
                outfile.write(line)

    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")

#Copies the metal-related lines for all sections except NONBON.
#Copies new atom type lines only for the ANGLE section.
def copy_metal_and_new_atom_type_lines(forcefield_file, output_file, metal_atom_types, new_atom_types):
    
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

            if current_section in section_lines:
                if current_section == "ANGLE":
                    # For ANGLE section, copy both metal and new atom type lines
                    if any(atom_type in stripped_line for atom_type in metal_atom_types.union(new_atom_types)):
                        section_lines[current_section].append(line)
                elif current_section in ["BOND", "DIHE", "IMPROPER"]:
                    # For BOND, DIHE, IMPROPER sections, only copy metal-related lines
                    if any(atom_type in stripped_line for atom_type in metal_atom_types):
                        section_lines[current_section].append(line)
                # For NONBON section, don't copy any lines

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
        print(f"Error inserting metal and new atom type lines into the correct sections: {str(e)}")


# Update the info for MASS section especially for metals
def update_mass_section(input_file, output_file, metal_atom_types, atom_info):
    try:
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

                    parts = line.split()
                    if len(parts) >= 3 and parts[0] in metal_atom_types:
                        atom_type = parts[0]
                        current_mass = float(parts[1])

                        # Find the corresponding atom_name for this atom_type
                        atom_name = next((info['name'] for info in atom_info.values() if info['type'] == atom_type), None)

                        if atom_name:
                            # Get the correct atomic mass from the periodictable library
                            element = getattr(periodictable, atom_name, None)
                            if element:
                                correct_mass = element.mass
                                if abs(current_mass - correct_mass) > 0.1:  # Allow for small differences due to rounding
                                    # Replace the mass with the correct value
                                    new_line = f"{atom_type} {correct_mass:.3f}{line[len(atom_type) + len(parts[1]) + 2:]}"
                                    outfile.write(new_line)
                                else:
                                    outfile.write(line)
                            else:
                                print(f"Warning: Element {atom_name} not found in periodictable library.")
                                outfile.write(line)
                        else:
                            print(f"Warning: No corresponding atom name found for atom type {atom_type}.")
                            outfile.write(line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)

    except Exception as e:
        print(f"Error updating MASS section: {str(e)}")

def main():
    """
    Main function that orchestrates the overall workflow.
    """
    atom_info = extract_atom_types_from_mol2('COMPLEX.mol2')

    # Read metal numbers
    metal_ids = read_metal_numbers('metal_number.dat')

    # Find atom types corresponding to metal IDs
    metal_atom_types = {atom_info[metal_id]['type'] for metal_id in metal_ids if metal_id in atom_info}

    # Read new atom types
    new_atom_types = read_new_atom_types('new_atomtype.dat')

    # Update FRCMOD file by removing lines related to metals and new atom types
    remove_metal_and_new_atom_type_lines('updated_COMPLEX_modified.frcmod', 'temp_COMPLEX_modified.frcmod', metal_atom_types, new_atom_types)

    # Update MASS section
    update_mass_section('temp_COMPLEX_modified.frcmod', 'updated_COMPLEX_modified2.frcmod', metal_atom_types, atom_info)

    # Copy metal-related and new atom type lines from forcefield.dat to updated_COMPLEX_modified2.frcmod in the correct sections
    # For BOND section, only copy metal-related lines
    copy_metal_and_new_atom_type_lines('forcefield2.dat', 'updated_COMPLEX_modified2.frcmod', metal_atom_types, new_atom_types)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()

