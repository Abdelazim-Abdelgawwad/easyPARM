import re

#READ the bond info from forcefield2.dat
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

#READ the angle info from forcefield2.dat
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

#READ the dihedral info from forcefield2.dat
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
#Remove any duplication or unnecessary files
def process_frcmod_file(reference_file, frcmod_file, output_file):
    # Read reference data
    reference_bonds = read_bond_data(reference_file)
    reference_angles = read_angle_data(reference_file)
    reference_dihes = read_dihe_data(reference_file)

    # Process and filter frcmod file
    bond_section = False
    angle_section = False
    dihe_section = False
    bond_pattern = re.compile(r"(\S+)\s*-\s*(\S+)")
    angle_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")
    dihe_pattern = re.compile(r"(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s*-\s*(\S+)")

    with open(frcmod_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            original_line = line
            line = line.split('#')[0].strip()  # Remove comments

            if line == "BOND":
                bond_section = True
                angle_section = False
                dihe_section = False
                outfile.write(original_line)
                continue
            elif line == "ANGLE":
                bond_section = False
                angle_section = True
                dihe_section = False
                outfile.write(original_line)
                continue
            elif line == "DIHE":
                bond_section = False
                angle_section = False
                dihe_section = True
                outfile.write(original_line)
                continue
            elif line == "IMPROPER":
                bond_section = False
                angle_section = False
                dihe_section = False
                outfile.write(original_line)
                continue

            if bond_section:
                match = bond_pattern.match(line)
                if match:
                    col1, col2 = match.groups()
                    bond = tuple(sorted([col1.strip(), col2.strip()]))
                    if bond in reference_bonds:
                        outfile.write(original_line)
                else:
                    outfile.write(original_line)
            elif angle_section:
                match = angle_pattern.match(line)
                if match:
                    col1, col2, col3 = match.groups()
                    angle1 = tuple(sorted([col1.strip(), col3.strip()]) + [col2.strip()])
                    angle2 = tuple(sorted([col3.strip(), col1.strip()]) + [col2.strip()])
                    if angle1 in reference_angles or angle2 in reference_angles:
                        outfile.write(original_line)
                else:
                    outfile.write(original_line)
            elif dihe_section:
                match = dihe_pattern.match(line)
                if match:
                    col1, col2, col3, col4 = match.groups()
                    dihe1 = tuple([col1.strip(), col2.strip(), col3.strip(), col4.strip()])
                    dihe2 = tuple([col4.strip(), col3.strip(), col2.strip(), col1.strip()])
                    if dihe1 in reference_dihes or dihe2 in reference_dihes:
                        outfile.write(original_line)
                else:
                    outfile.write(original_line)
            else:
                outfile.write(original_line)

# Call the function 
reference_file = 'forcefield2.dat'
frcmod_file = 'updated_updated_COMPLEX_modified2.frcmod'
output_file = 'filtered_COMPLEX_modified2.frcmod'
process_frcmod_file(reference_file, frcmod_file, output_file)
