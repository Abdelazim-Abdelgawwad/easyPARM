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


import re

def read_bond_data(file_path):
    bond_data = set()
    bond_section = False
    pattern = re.compile(r"(\S+)\s*-\s*(\S+)")
    with open(file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip()  # Remove comments
            if line == "BOND":
                bond_section = True
                continue
            elif line == "ANGLE":
                break
            if bond_section:
                match = pattern.match(line)
                if match:
                    col1, col2 = match.groups()
                    bond = tuple(sorted([col1.strip(), col2.strip()]))
                    bond_data.add(bond)
    return bond_data

def read_angle_data(file_path):
    angle_data = set()
    angle_section = False
    pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    with open(file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip()  # Remove comments
            if line == "ANGLE":
                angle_section = True
                continue
            elif line == "DIHE":
                break
            if angle_section:
                match = pattern.match(line)
                if match:
                    col1, col2, col3 = match.groups()
                    angle = tuple(sorted([col1.strip(), col3.strip()]) + [col2.strip()])
                    angle_data.add(angle)
    return angle_data

def read_dihe_data(file_path):
    dihe_data = set()
    dihe_section = False
    pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    with open(file_path, 'r') as file:
        for line in file:
            line = line.split('#')[0].strip()  # Remove comments
            if line == "DIHE":
                dihe_section = True
                continue
            elif line == "IMPROPER":
                break
            if dihe_section:
                match = pattern.match(line)
                if match:
                    col1, col2, col3, col4 = match.groups()
                    dihe = tuple([col1.strip(), col2.strip(), col3.strip(), col4.strip()])
                    dihe_reversed = tuple([col4.strip(), col3.strip(), col2.strip(), col1.strip()])
                    dihe_data.add(dihe)
                    dihe_data.add(dihe_reversed)
    return dihe_data

#Clean MASS section by consolidating entries.
def clean_mass_section(lines):
    mass_entries = {}
    cleaned_lines = []
    in_mass_section = False
    pattern = re.compile(r"(\S+)\s+(\S+)\s+(\S+)")

    # First pass: collect all unique mass entries
    for line in lines:
        stripped_line = line.strip()
        if stripped_line == "MASS":
            in_mass_section = True
            continue
        elif stripped_line in ["BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
            in_mass_section = False
            continue
            
        if in_mass_section:
            match = pattern.match(stripped_line)
            if match:
                atom_type = match.group(1).strip()
                mass = float(match.group(2))
                polarizability = float(match.group(3))
                
                # Store the entry with highest polarizability for each atom type
                if atom_type not in mass_entries or polarizability > mass_entries[atom_type][1]:
                    mass_entries[atom_type] = (mass, polarizability, line)

    # Second pass: reconstruct the file with unique entries
    in_mass_section = False
    for line in lines:
        stripped_line = line.strip()
        if stripped_line == "MASS":
            in_mass_section = True
            cleaned_lines.append("MASS\n")  # Add MASS header
            # Add all unique mass entries after MASS header
            for _, (_, _, original_line) in sorted(mass_entries.items()):
                cleaned_lines.append(original_line)
            cleaned_lines.append("\n")  # Add blank line after MASS section
            continue
        elif stripped_line in ["BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
            in_mass_section = False
            cleaned_lines.append(line)
            continue
            
        if not in_mass_section:
            cleaned_lines.append(line)

    return cleaned_lines

def process_and_clean_frcmod_file(reference_file, frcmod_file, output_file):
    # Read reference data
    reference_bonds = read_bond_data(reference_file)
    reference_angles = read_angle_data(reference_file)
    reference_dihes = read_dihe_data(reference_file)

    # First read all lines and clean MASS section
    with open(frcmod_file, 'r') as infile:
        lines = infile.readlines()
    
    # Clean MASS section first
    lines = clean_mass_section(lines)

    # Process and filter other sections
    bond_section = False
    angle_section = False
    dihe_section = False
    nonbon_section = False
    
    # Patterns for sections
    bond_pattern = re.compile(r"(\S+)\s*-\s*(\S+)")
    angle_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    dihe_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    nonbon_pattern = re.compile(r"(\S+)\s+(\S+)\s+(\S+)")

    # Store unique NONBON entries
    nonbon_entries = {}
    output_lines = []
    last_section = None

    for line in lines:
        original_line = line
        stripped_line = line.split('#')[0].strip()

        # Section detection
        if stripped_line == "BOND":
            if last_section != "MASS":  # Only add blank line if not following MASS section
                output_lines.append("\n")
            bond_section = True
            angle_section = dihe_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "BOND"
            continue
        elif stripped_line == "ANGLE":
            output_lines.append("\n")
            angle_section = True
            bond_section = dihe_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "ANGLE"
            continue
        elif stripped_line == "DIHE":
            output_lines.append("\n")
            dihe_section = True
            bond_section = angle_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "DIHE"
            continue
        elif stripped_line == "NONBON":
            output_lines.append("\n")
            nonbon_section = True
            bond_section = angle_section = dihe_section = False
            output_lines.append(original_line)
            last_section = "NONBON"
            continue
        elif stripped_line == "IMPROPER":
            output_lines.append("\n")
            bond_section = angle_section = dihe_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "IMPROPER"
            continue

        # Process each section
        if bond_section:
            match = bond_pattern.match(stripped_line)
            if match:
                col1, col2 = match.groups()
                bond = tuple(sorted([col1.strip(), col2.strip()]))
                if bond in reference_bonds:
                    output_lines.append(original_line)
            else:
                output_lines.append(original_line)

        elif angle_section:
            match = angle_pattern.match(stripped_line)
            if match:
                col1, col2, col3 = match.groups()
                angle1 = tuple(sorted([col1.strip(), col3.strip()]) + [col2.strip()])
                angle2 = tuple(sorted([col3.strip(), col1.strip()]) + [col2.strip()])
                if angle1 in reference_angles or angle2 in reference_angles:
                    output_lines.append(original_line)
            else:
                output_lines.append(original_line)

        elif dihe_section:
            match = dihe_pattern.match(stripped_line)
            if match:
                col1, col2, col3, col4 = match.groups()
                dihe1 = tuple([col1.strip(), col2.strip(), col3.strip(), col4.strip()])
                dihe2 = tuple([col4.strip(), col3.strip(), col2.strip(), col1.strip()])
                if dihe1 in reference_dihes or dihe2 in reference_dihes:
                    output_lines.append(original_line)
            else:
                output_lines.append(original_line)

        elif nonbon_section:
            match = nonbon_pattern.match(stripped_line)
            if match:
                atom_type, param1, param2 = match.groups()
                key = atom_type.strip()
                if key not in nonbon_entries:
                    nonbon_entries[key] = original_line
                    output_lines.append(original_line)
            else:
                output_lines.append(original_line)

        else:
            output_lines.append(original_line)

    # Write the cleaned and filtered output
    with open(output_file, 'w') as outfile:
        outfile.writelines(output_lines)

# input
reference_file = 'forcefield2.dat'
frcmod_file = 'updated_updated_COMPLEX_modified2.frcmod'
output_file = 'filtered_COMPLEX_modified2.frcmod'
process_and_clean_frcmod_file(reference_file, frcmod_file, output_file)
