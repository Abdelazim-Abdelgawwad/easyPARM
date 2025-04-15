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
#                              |  $$$$$$/              Ver. 3.25 - 14 April 2025                                  #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################



from collections import defaultdict

# Function to process the BOND and ANGLE sections
def process_bond_angle_data(data):
    # Use defaultdict to automatically initialize new keys with default values
    bond_angle_data = defaultdict(lambda: [0.0, 0.0, 0])
    for line in data:
        if line.strip():  # Skip empty lines
            fields = line.split()
            atom_types = fields[0].split('-')
            key = '-'.join(atom_types)
            # Accumulate values and count for each key
            bond_angle_data[key][0] += float(fields[1])  # Sum of first values
            bond_angle_data[key][1] += float(fields[2])  # Sum of second values
            bond_angle_data[key][2] += 1  # Count of entries
    return bond_angle_data

# Read the input file
with open('forcefield.dat', 'r') as file:
    content = file.readlines()

# Find the start and end of the BOND, ANGLE, and DIHE sections
bond_start = None
bond_end = None
angle_start = None
angle_end = None
dihe_start = None
dihe_end = None

# Iterate through the file to find section boundaries
for i, line in enumerate(content):
    if line.strip() == 'BOND':
        bond_start = i + 1
    elif line.strip() == 'ANGLE':
        bond_end = i
        angle_start = i + 1
    elif line.strip() == 'DIHE':
        angle_end = i
        dihe_start = i
        # Find the end of DIHE section
        for j, dihe_line in enumerate(content[i:]):
            if dihe_line.strip() == 'IMPROPER':
                dihe_end = i + j
                break

# Process the BOND and ANGLE sections
bond_data = process_bond_angle_data(content[bond_start:bond_end])
angle_data = process_bond_angle_data(content[angle_start:angle_end])

# Write the updated data to a new file
with open('forcefield2.dat', 'w') as file:
    # Write content before BOND section
    for line in content[:bond_start]:
        file.write(line)
    
    # Write processed BOND data
    for key, value in bond_data.items():
        file.write(f"{key} {value[0] / value[2]:.2f} {value[1] / value[2]:.3f}\n")
    file.write("\n")  # Add a blank line after the BOND section
    
    # Write content between BOND and ANGLE sections
    for line in content[bond_end:angle_start]:
        file.write(line)
    
    # Write processed ANGLE data
    for key, value in angle_data.items():
        file.write(f"{key} {value[0] / value[2]:.2f} {value[1] / value[2]:.3f}\n")
    file.write("\n")  # Add a blank line after the ANGLE section
    
    # Write DIHE section
    for line in content[dihe_start:dihe_end]:
        file.write(line)
    file.write("\n")  # Add a blank line after the DIHE section
    
    # Write remaining content
    for line in content[dihe_end:]:
        file.write(line)
