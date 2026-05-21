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
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
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

#Check if all three bond pairs from an improper exist in the reference bonds.
def check_improper_bonds_exist(col1, col2, col3, col4, reference_bonds):
    bond1 = tuple(sorted([col1.strip(), col2.strip()]))
    bond2 = tuple(sorted([col2.strip(), col3.strip()]))
    bond3 = tuple(sorted([col3.strip(), col4.strip()]))
    return bond1 in reference_bonds and bond2 in reference_bonds and bond3 in reference_bonds

#Clean MASS section.
def clean_mass_section(lines):
    mass_entries = {}
    cleaned_lines = []
    in_mass_section = False
    pattern = re.compile(r"(\S+)\s+(\S+)\s+(\S+)")

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
                if atom_type not in mass_entries or polarizability > mass_entries[atom_type][1]:
                    mass_entries[atom_type] = (mass, polarizability, line)

    in_mass_section = False
    for line in lines:
        stripped_line = line.strip()
        if stripped_line == "MASS":
            in_mass_section = True
            cleaned_lines.append("MASS\n")
            for _, (_, _, original_line) in sorted(mass_entries.items()):
                cleaned_lines.append(original_line)
            cleaned_lines.append("\n")
            continue
        elif stripped_line in ["BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
            in_mass_section = False
            cleaned_lines.append(line)
            continue
            
        if not in_mass_section:
            cleaned_lines.append(line)

    return cleaned_lines

def process_and_clean_frcmod_file(reference_file, frcmod_file, output_file):
    reference_bonds = read_bond_data(reference_file)
    reference_angles = read_angle_data(reference_file)
    reference_dihes = read_dihe_data(reference_file)

    with open(frcmod_file, 'r') as infile:
        lines = infile.readlines()
    
    lines = clean_mass_section(lines)

    bond_section = False
    angle_section = False
    dihe_section = False
    improper_section = False
    nonbon_section = False
    
    bond_pattern = re.compile(r"(\S+)\s*-\s*(\S+)")
    angle_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    dihe_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    improper_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    nonbon_pattern = re.compile(r"(\S+)\s+(\S+)\s+(\S+)")

    nonbon_entries = {}
    output_lines = []
    last_section = None

    for line in lines:
        original_line = line
        stripped_line = line.split('#')[0].strip()

        if stripped_line == "BOND":
            if last_section != "MASS":
                output_lines.append("\n")
            bond_section = True
            angle_section = dihe_section = improper_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "BOND"
            continue
        elif stripped_line == "ANGLE":
            output_lines.append("\n")
            angle_section = True
            bond_section = dihe_section = improper_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "ANGLE"
            continue
        elif stripped_line == "DIHE":
            output_lines.append("\n")
            dihe_section = True
            bond_section = angle_section = improper_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "DIHE"
            continue
        elif stripped_line == "IMPROPER":
            output_lines.append("\n")
            improper_section = True
            bond_section = angle_section = dihe_section = nonbon_section = False
            output_lines.append(original_line)
            last_section = "IMPROPER"
            continue
        elif stripped_line == "NONBON":
            output_lines.append("\n")
            nonbon_section = True
            bond_section = angle_section = dihe_section = improper_section = False
            output_lines.append(original_line)
            last_section = "NONBON"
            continue

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

        elif improper_section:
            match = improper_pattern.match(stripped_line)
            if match:
                col1, col2, col3, col4 = match.groups()
                if check_improper_bonds_exist(col1, col2, col3, col4, reference_bonds):
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

    with open(output_file, 'w') as outfile:
        outfile.writelines(output_lines)


#Fix atom-type padding so every type is exactly 2 characters.
#Return a 2-character atom-type string, padding with a trailing space if needed.
def pad_atom_type(token: str) -> str:
    t = token.strip()
    if len(t) == 1 and t.isalpha():
        return t + ' '
    return t


#Read *input_file* (a filtered .frcmod), rewrite every BOND / ANGLE /
def fix_atom_type_spacing(input_file: str, output_file: str) -> None:

    # How many '-'-separated tokens each section has
    SECTION_ATOM_COUNT = {
        "BOND":     2,
        "ANGLE":    3,
        "DIHE":     4,
        "IMPROPER": 4,
    }
    SECTION_KEYWORDS = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}

    # Regex: captures the leading type-string (tokens separated by '-')
    def build_type_re(n_atoms: int) -> re.Pattern:
        token = r"(\S+)"
        type_part = r"\s*-\s*".join([token] * n_atoms)
        return re.compile(r"^(" + type_part + r")(.*)", re.DOTALL)

    current_section = None
    output_lines = []

    with open(input_file, 'r') as fh:
        for line in fh:
            # Split off inline comment to check the keyword, but keep original
            keyword = line.split('#')[0].strip()

            #  section header 
            if keyword in SECTION_KEYWORDS:
                current_section = keyword
                output_lines.append(line)
                continue

            # blank / comment-only / non-data lines 
            if not keyword or current_section not in SECTION_ATOM_COUNT:
                output_lines.append(line)
                continue

            #  data line inside BOND / ANGLE / DIHE / IMPROPER
            n = SECTION_ATOM_COUNT[current_section]

            # Build a regex that captures each atom-type token individually
            # e.g. for BOND: ^(\S+)\s*-\s*(\S+)(.*)
            token_groups = r"\s*-\s*".join([r"(\S+)"] * n)
            data_re = re.compile(r"^" + token_groups + r"(.*)", re.DOTALL)

            match = data_re.match(line)
            if not match:
                # Doesn't look like a data line (e.g. a stray comment) → keep as-is
                output_lines.append(line)
                continue

            groups = match.groups()
            atom_types = groups[:n]       # the n atom-type tokens
            remainder  = groups[n]        # everything after the last token

            # Pad each single-letter atom type
            padded = [pad_atom_type(t) for t in atom_types]

            # Reconstruct: "XX-XX" with no extra spaces around the dashes
            type_string = "-".join(padded)

            new_line = type_string + remainder
            output_lines.append(new_line)

    with open(output_file, 'w') as fh:
        fh.writelines(output_lines)


# Main pipeline

# Step 1 – filter frcmod against the reference force-field
reference_file  = 'forcefield2.dat'
frcmod_file     = 'updated_updated_COMPLEX_modified2.frcmod'
filtered_file   = 'filtered_COMPLEX_modified2.frcmod'

process_and_clean_frcmod_file(reference_file, frcmod_file, filtered_file)

# Step 2 – fix atom-type padding in the filtered file
spaced_file = 'spaced_filtered_COMPLEX_modified2.frcmod'
fix_atom_type_spacing(filtered_file, spaced_file)

