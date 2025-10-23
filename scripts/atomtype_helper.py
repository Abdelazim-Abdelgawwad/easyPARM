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
#                              |  $$$$$$/              Ver. 4.15 - 17 October 2025                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

#Extract bond data from a mol2 file.
def read_mol2_bond_data(mol2_filename):
    bond_data = []
    sequence_number = 0
    
    with open(mol2_filename, 'r') as mol2_file:
        is_bond_section = False
        for line in mol2_file:
            line = line.strip()
            # Identify the bond section
            if line == "@<TRIPOS>BOND":
                is_bond_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                is_bond_section = False
                continue
            
            # Process lines in the bond section
            if is_bond_section and line:
                parts = line.split()
                if len(parts) >= 4:  # Ensure we have enough parts
                    try:
                        seq_num = int(parts[0])
                        atom1 = int(parts[1])
                        atom2 = int(parts[2])
                        
                        # Handle both integer and 'ar' cases for bond type
                        if parts[3] == 'ar':
                            bond_type = 'ar'
                        else:
                            bond_type = int(parts[3])
                            
                        bond_data.append((seq_num, atom1, atom2, bond_type))
                        sequence_number = seq_num  # Keep track of the last sequence number
                    except ValueError as e:
                        print(f"Error processing line: {line} - {e}")
    
    return bond_data, sequence_number

#Read bond information from distance_type.dat file.
def read_distance_type_data(dat_filename):
    dat_bonds = {}
    
    with open(dat_filename, 'r') as dat_file:
        for line in dat_file:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        bond_type = int(parts[2])
                        
                        # Store both directions of the bond for easier lookup
                        dat_bonds[(atom1, atom2)] = bond_type
                        dat_bonds[(atom2, atom1)] = bond_type
                    except ValueError as e:
                        print(f"Error processing line: {line} - {e}")
    
    return dat_bonds

#Update bond information in mol2 file with data from distance_type.dat file.
def update_mol2_bonds(mol2_filename, dat_filename, output_filename):
    # Read bond data from both files
    mol2_bonds, last_seq = read_mol2_bond_data(mol2_filename)
    dat_bonds = read_distance_type_data(dat_filename)
    
    print(f"Found {len(mol2_bonds)} bonds in mol2 file")
    print(f"Found {len(dat_bonds)//2} unique bonds in dat file (stored bidirectionally)")
    
    # Create a dictionary to track updated bonds
    updated_bonds = {}
    updated_count = 0
    
    # Read the original mol2 file
    with open(mol2_filename, 'r') as mol2_file:
        mol2_content = mol2_file.readlines()
    
    # Process the mol2 file
    with open(output_filename, 'w') as output_file:
        is_bond_section = False
        for line in mol2_content:
            original_line = line
            line_stripped = line.strip()
            
            # Track if we're in the bond section
            if line_stripped == "@<TRIPOS>BOND":
                is_bond_section = True
                output_file.write(original_line)
                continue
            if line_stripped.startswith("@<TRIPOS>") and line_stripped != "@<TRIPOS>BOND":
                is_bond_section = False
                output_file.write(original_line)
                continue
            
            # Update bond information
            if is_bond_section and line_stripped:
                parts = line_stripped.split()
                if len(parts) >= 4:
                    try:
                        seq_num = int(parts[0])
                        atom1 = int(parts[1])
                        atom2 = int(parts[2])
                        
                        # Check if this bond exists in the dat file
                        if (atom1, atom2) in dat_bonds:
                            new_bond_type = dat_bonds[(atom1, atom2)]
                            # Create updated line with proper formatting to match mol2 format
                            # Preserve the original spacing/format as much as possible
                            updated_line = f"     {seq_num:<6d}{atom1:<6d}{atom2:<6d}{new_bond_type:<3d}\n"
                            output_file.write(updated_line)
                            updated_bonds[(atom1, atom2)] = True
                            updated_count += 1
                        else:
                            # Keep original bond information
                            output_file.write(original_line)
                    except ValueError as e:
                        # If there's an error parsing, keep the original line
                        print(f"Error processing line: {line_stripped} - {e}")
                        output_file.write(original_line)
                else:
                    # Write the original line if it doesn't match expected format
                    output_file.write(original_line)
            else:
                # Write non-bond section lines as is
                output_file.write(original_line)
    
    return updated_count

def main():
    import sys
    
    if len(sys.argv) != 4:
        print("Usage: python3 atomtype_helper.py input.mol2 input.dat output.mol2")
        sys.exit(1)
        
    mol2_filename = sys.argv[1]
    dat_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    print(f"Reading mol2 file: {mol2_filename}")
    print(f"Reading dat file: {dat_filename}")
    print(f"Output will be written to: {output_filename}")
    
    updated_count = update_mol2_bonds(mol2_filename, dat_filename, output_filename)
    print(f"Process complete. Updated {updated_count} bonds.")

if __name__ == "__main__":
    main()
