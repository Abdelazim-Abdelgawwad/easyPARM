import re

# Read the atom types from the file
def process_frcmod_file(file_path, atomtype_file):
    with open(atomtype_file, 'r') as at_file:
        atom_types = set(line.strip() for line in at_file)

    # Read the input file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize section flags and data dictionaries
    bond_section = False
    angle_section = False
    dihe_section = False
    data_bond = {}
    data_angle = {}
    data_dihe = {}

    # Regex patterns for BOND, ANGLE, and DIHE sections
    pattern_bond = re.compile(r"(\S+)-(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
    pattern_angle = re.compile(r"(\S+)-(\S+)-(\S+)\s+(\d+\.\d+)\s+(\d+\.\d+)")
    pattern_dihe = re.compile(r"(\S+)-(\S+)-(\S+)-(\S+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")

    # Process the lines to extract data between BOND, ANGLE, and DIHE sections
    for line in lines:
        line = line.strip()
        
        if line == "BOND":
            bond_section = True
            angle_section = False
            dihe_section = False
            continue
        elif line == "ANGLE":
            bond_section = False
            angle_section = True
            dihe_section = False
            continue
        elif line == "DIHE":
            bond_section = False
            angle_section = False
            dihe_section = True
            continue

        # Process BOND section
        if bond_section:
            match = pattern_bond.match(line)
            if match:
                col1, col2, col3, col4 = match.groups()
                col3, col4 = float(col3), float(col4)

                key = f"{col1}-{col2}"
                if key not in data_bond:
                    data_bond[key] = {'sum_col3': 0.0, 'sum_col4': 0.0, 'count': 0}
                data_bond[key]['sum_col3'] += col3
                data_bond[key]['sum_col4'] += col4
                data_bond[key]['count'] += 1

        # Process ANGLE section
        elif angle_section:
            match = pattern_angle.match(line)
            if match:
                col1, col2, col3, col4, col5 = match.groups()
                col4, col5 = float(col4), float(col5)

                key = f"{col1}-{col2}-{col3}"
                if key not in data_angle:
                    data_angle[key] = {'sum_col4': 0.0, 'sum_col5': 0.0, 'count': 0}
                data_angle[key]['sum_col4'] += col4
                data_angle[key]['sum_col5'] += col5
                data_angle[key]['count'] += 1

        # Process DIHE section
        elif dihe_section:
            match = pattern_dihe.match(line)
            if match:
                col1, col2, col3, col4, col5, col6, col7, col8 = match.groups()
                col5 = int(col5)
                col6, col7, col8 = float(col6), float(col7), float(col8)

                key = f"{col1}-{col2}-{col3}-{col4}"
                if key not in data_dihe:
                    data_dihe[key] = {'sum_col5': 0, 'sum_col6': 0.0, 'sum_col7': 0.0, 'sum_col8': 0.0, 'count': 0}
                data_dihe[key]['sum_col5'] += col5
                data_dihe[key]['sum_col6'] += col6
                data_dihe[key]['sum_col7'] += col7
                data_dihe[key]['sum_col8'] += col8
                data_dihe[key]['count'] += 1

    # Write the updated content to a new file
    output_file_path = 'updated_' + file_path
    with open(output_file_path, 'w') as output_file:
        bond_section = False
        angle_section = False
        dihe_section = False
        for line in lines:
            line = line.strip()
            # Set section flags based on section headers
            if line == "BOND":
                bond_section = True
                angle_section = False
                dihe_section = False
                output_file.write(line + "\n")
                continue
            elif line == "ANGLE":
                bond_section = False
                angle_section = True
                dihe_section = False
                output_file.write(line + "\n")
                continue
            elif line == "DIHE":
                bond_section = False
                angle_section = False
                dihe_section = True
                output_file.write(line + "\n")
                continue

            # Process and write BOND section
            if bond_section:
                match = pattern_bond.match(line)
                if match:
                    col1, col2, _, _ = match.groups()
                    key = f"{col1}-{col2}"
                    if key in data_bond:
                        sum_col3 = data_bond[key]['sum_col3']
                        sum_col4 = data_bond[key]['sum_col4']
                        count = data_bond[key]['count']
                        avg_col3 = sum_col3 / count
                        avg_col4 = sum_col4 / count
                        output_file.write(f"{col1}-{col2} {avg_col3:.2f} {avg_col4:.3f}\n")
                        data_bond.pop(key)  # Remove the key after processing
                else:
                    output_file.write(line + "\n")

            # Process and write ANGLE section
            elif angle_section:
                match = pattern_angle.match(line)
                if match:
                    col1, col2, col3, _, _ = match.groups()
                    key = f"{col1}-{col2}-{col3}"
                    if key in data_angle:
                        sum_col4 = data_angle[key]['sum_col4']
                        sum_col5 = data_angle[key]['sum_col5']
                        count = data_angle[key]['count']
                        avg_col4 = sum_col4 / count
                        avg_col5 = sum_col5 / count
                        output_file.write(f"{col1}-{col2}-{col3} {avg_col4:.3f} {avg_col5:.3f}\n")
                        data_angle.pop(key)  # Remove the key after processing
                else:
                    output_file.write(line + "\n")

            # Process and write DIHE section
            elif dihe_section:
                match = pattern_dihe.match(line)
                if match:
                    col1, col2, col3, col4, _, _, _, _ = match.groups()
                    key = f"{col1}-{col2}-{col3}-{col4}"
                    if key in data_dihe:
                        sum_col5 = data_dihe[key]['sum_col5']
                        sum_col6 = data_dihe[key]['sum_col6']
                        sum_col7 = data_dihe[key]['sum_col7']
                        sum_col8 = data_dihe[key]['sum_col8']
                        count = data_dihe[key]['count']
                        avg_col5 = round(sum_col5 / count)
                        avg_col6 = sum_col6 / count
                        avg_col7 = sum_col7 / count
                        avg_col8 = sum_col8 / count
                        output_file.write(f"{col1}-{col2}-{col3}-{col4}   {avg_col5:d}    {avg_col6:.3f}       {avg_col7:.4f}           {avg_col8:.3f}\n")
                        data_dihe.pop(key)  # Remove the key after processing
                else:
                    output_file.write(line + "\n")

            # Write lines that don't belong to any section as it is
            else:
                output_file.write(line + "\n")

# Call the function with the file paths
process_frcmod_file('updated_COMPLEX_modified2.frcmod', 'new_atomtype.dat')
