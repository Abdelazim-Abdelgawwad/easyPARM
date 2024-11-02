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


from collections import defaultdict

#Extract atom types from a MOL2 file.
def extract_atom_types_from_mol2(file_path):
        
    with open(file_path, 'r') as file:
        lines = file.readlines()
    atom_types = {}
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
                atom_type = parts[5]
                atom_types[atom_id] = atom_type
    return atom_types

#Process bond, angle, and dihedral information from an input file and write to an output file.
def process_bond_angle_dihedral_file(atom_types, input_file, output_file):
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        section = None
        
        # Write header
        outfile.write("Remark line goes here\n")
        outfile.write("MASS\n")
        
        for line in infile:
            line = line.strip()
            
            # Determine the current section (BOND, ANGLE, or DIHE)
            if line.startswith("Bond Information:"):
                section = "BOND"
                outfile.write("\n")
                outfile.write("BOND\n")
                continue
            elif line.startswith("Angle Information:"):
                section = "ANGLE"
                outfile.write("\n")
                outfile.write("ANGLE\n")
                continue
            elif line.startswith("Dihedral Information:"):
                section = "DIHE"
                outfile.write("\n")
                outfile.write("DIHE\n")
                continue
            
            parts = line.split()
            
            # Process BOND information
            if section == "BOND" and len(parts) == 4:
                atom1_id, atom2_id, distance, force = int(parts[0]), int(parts[1]), parts[2], parts[3]
                atom1_type = atom_types.get(atom1_id, "Unknown")
                atom2_type = atom_types.get(atom2_id, "Unknown")
                outfile.write(f"{atom1_type}-{atom2_type}    {force}    {distance}\n")
            
            # Process ANGLE information
            elif section == "ANGLE" and len(parts) == 5:
                atom1_id, atom2_id, atom3_id, angle, force = int(parts[0]), int(parts[1]), int(parts[2]), parts[3], parts[4]
                atom1_type = atom_types.get(atom1_id, "Unknown")
                atom2_type = atom_types.get(atom2_id, "Unknown")
                atom3_type = atom_types.get(atom3_id, "Unknown")
                outfile.write(f"{atom1_type}-{atom2_type}-{atom3_type}     {force}     {angle}\n")
            
            # Process DIHE (dihedral) information
            elif section == "DIHE" and len(parts) == 5:
                atom1_id, atom2_id, atom3_id, atom4_id, dihedral = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), parts[4]
                atom1_type = atom_types.get(atom1_id, "Unknown")
                atom2_type = atom_types.get(atom2_id, "Unknown")
                atom3_type = atom_types.get(atom3_id, "Unknown")
                atom4_type = atom_types.get(atom4_id, "Unknown")
                outfile.write(f"{atom1_type}-{atom2_type}-{atom3_type}-{atom4_type}   1    0.000       0.0000           1.000\n")
        
        # Write footer
        outfile.write("\n")
        outfile.write("IMPROPER\n")
        outfile.write("\n")
        outfile.write("NONBON\n")
        outfile.write("\n")
            
# Inputs
mol2_file_path = 'NEW_COMPLEX.mol2'
input_file_path = 'bond_angle_dihedral_data.dat'
# Output
output_file_path = 'forcefield.dat'

# Extract atom types from the MOL2 file
atom_types = extract_atom_types_from_mol2(mol2_file_path)

# Process the bond, angle, and dihedral file and save the results
process_bond_angle_dihedral_file(atom_types, input_file_path, output_file_path)
