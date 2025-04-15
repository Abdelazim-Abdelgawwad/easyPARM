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


import re
import os

#Extract atom types from mol2 file and save them to Reference_atom_type.dat
def extract_atom_types_from_mol2(file_path):
    try:
        # First read line_number.dat to get the replacements
        replacements = {}
        with open('line_number.dat', 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    line_num = int(parts[0])
                    old_atom = parts[1]
                    # Remove any numbers from the old atom name
                    base_atom = ''.join(char for char in old_atom if not char.isdigit())
                    replacements[line_num] = base_atom

        # Read the mol2 file
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        atom_info = {}
        is_atom_section = False
        atom_types = []
        
        # Process mol2 file and collect atom types
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
                    atom_types.append(atom_type)
        
        # Apply replacements based on line_number.dat
        for line_num, base_atom in replacements.items():
            if 1 <= line_num <= len(atom_types):
                atom_types[line_num - 1] = base_atom
        
        # Save modified atom types to Reference_atom_type.dat
        with open('Reference_atom_type.dat', 'w') as f:
            for atom_type in atom_types:
                f.write(f"{atom_type}\n")
                
        return atom_info
    except Exception as e:
        print(f"Error processing files: {str(e)}")
        return {}

#Extract bond section data from mol2 file and return both lines and bond types
def extract_bond_section(mol2_filename):
    bond_section = []
    bond_types = []
    with open(mol2_filename, 'r') as mol2_file:
        is_bond_section = False
        for line in mol2_file:
            if line.startswith("@<TRIPOS>BOND"):
                is_bond_section = True
                bond_section.append(line)
                continue
            if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                break
            if is_bond_section:
                parts = line.split()
                if len(parts) >= 4:
                    bond_types.append(parts[3].lower())  # Store bond type
                bond_section.append(line)
    return bond_section, bond_types

#Replace atom types in COMPLEX.mol2 with data from ONE.mol2, and only replace bonds if no 'ar' bonds are present
def replace_atom_types_and_bonds_from_mol2(complex_file='COMPLEX.mol2', reference_file='ONE.mol2'):
    try:
        # First check if COMPLEX.mol2 has any 'ar' bonds
        complex_bonds, complex_bond_types = extract_bond_section(complex_file)
        has_aromatic = any(bond_type == 'ar' for bond_type in complex_bond_types)
        
        # Read Reference_atom_type.dat
        with open('Reference_atom_type.dat', 'r') as f:
            reference_types = [line.strip() for line in f.readlines()]

        # Get bond section from reference file only if we'll need it
        reference_bonds = None
        if not has_aromatic:
            reference_bonds, _ = extract_bond_section(reference_file)
        
        # Read COMPLEX.mol2
        with open(complex_file, 'r') as f:
            lines = f.readlines()

        new_lines = []
        is_atom_section = False
        is_bond_section = False
        atom_index = 0

        # Process COMPLEX.mol2
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                is_bond_section = False
                new_lines.append(line)
                continue
            elif line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                is_bond_section = True
                new_lines.append(line)
                if has_aromatic:
                    # Keep original bonds if aromatic
                    new_lines.extend(complex_bonds[1:])  # Skip the @<TRIPOS>BOND header
                else:
                    # Replace with reference bonds
                    new_lines.extend(reference_bonds[1:])  # Skip the @<TRIPOS>BOND header
                continue
            elif line.startswith("@<TRIPOS>SUBSTRUCTURE"):
                is_atom_section = False
                is_bond_section = False
                new_lines.append(line)
                continue
            
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6 and atom_index < len(reference_types):
                    # Preserve the original line formatting by replacing just the atom type
                    prefix = line[:line.rfind(parts[5])]
                    suffix = line[line.rfind(parts[5]) + len(parts[5]):]
                    new_line = prefix + reference_types[atom_index] + suffix
                    new_lines.append(new_line)
                    atom_index += 1
                else:
                    new_lines.append(line)
            elif not is_bond_section:  # Only add non-bond section lines
                new_lines.append(line)

        # Write modified content to new file
        output_file = 'COMPLEX_modified.mol2'
        with open(output_file, 'w') as f:
            f.writelines(new_lines)
        
        return True

    except Exception as e:
        print(f"Error replacing atom types and bonds: {str(e)}")
        return False

#Process PDB file to replace N-1 metal atoms with S, where N is total number of metals. If only one metal is present, it will be replaced.
def process_pdb(input_file='input.pdb'):
    # Dictionary of metals and their atomic numbers
    metals = {
            'al': 13, 'ar': 18, 'cr': 24, 'sc': 21, 'ti': 22, 'v': 23, 'mn': 25, 'b':5,
        'fe': 26, 'co': 27, 'ni': 28, 'cu': 29, 'zn': 30, 'ga': 31, 'se': 34,
        'kr': 36, 'rb': 37, 'sr': 38, 'y': 39, 'zr': 40, 'mo': 42,
        'tc': 43, 'ru': 44, 'rh': 45, 'pd': 46, 'ag': 47, 'cd': 48, 'in': 49,
        'sn': 50, 'xe': 54, 'cs': 55, 'ba': 56, 'la': 57, 'ce': 58, 'pr': 59,
        'nd': 60, 'pm': 61, 'sm': 62, 'eu': 63, 'gd': 64, 'tb': 65, 'dy': 66,
        'ho': 67, 'er': 68, 'tm': 69, 'yb': 70, 'lu': 71, 'hf': 72, 'ta': 73,
        'w': 74, 're': 75, 'os': 76, 'ir': 77, 'pt': 78, 'au': 79, 'hg': 80,
        'tl': 81, 'pb': 82, 'bi': 83, 'po': 84, 'at': 85, 'rn': 86, 'fr': 87,
        'ra': 88, 'ac': 89, 'th': 90, 'pa': 91, 'u': 92, 'np': 93, 'pu': 94,
        'am': 95, 'cm': 96, 'bk': 97, 'cf': 98, 'es': 99, 'fm': 100, 'md': 101,
        'no': 102, 'lr': 103, 'si': 14, 'ge': 32, 'as': 33, 'sb': 51, 'te': 52
    }

    # First pass: count total metals
    metal_atoms = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        
    for line_num, line in enumerate(lines, 1):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_name = line[12:16].strip()
            match = re.match(r'([A-Za-z]+)(\d*)', atom_name)
            if match:
                base_atom, _ = match.groups()
                base_atom_lower = base_atom.lower()
                if base_atom_lower in metals:
                    metal_atoms.append((line_num, line, atom_name))

    total_metals = len(metal_atoms)
    if total_metals == 0:
        with open('no_metal.dat', 'w') as f:
            f.write("No Metal\n")
        return

    if total_metals == 1:
        metals_to_replace = 1
    else:
        metals_to_replace = total_metals - 1
    
    if total_metals >= 5:
        with open('limited_data.dat', 'w') as f:
            for _, _, atom_name in metal_atoms:
                f.write(f"{atom_name}\n")

    metal_count = {}
    line_mappings = []
    new_lines = lines.copy()
    
    for i, (line_num, line, atom_name) in enumerate(metal_atoms):
        if i >= metals_to_replace:
            break
            
        match = re.match(r'([A-Za-z]+)(\d*)', atom_name)
        if match:
            base_atom, number = match.groups()
            number = number if number else '1'
            old_atom = f"{base_atom}{number}"
            
            if old_atom not in metal_count:
                metal_count[old_atom] = len(metal_count) + 1
            
            new_atom = f"S{metal_count[old_atom]}"
            line_mappings.append(f"{line_num} {old_atom} {new_atom}")
            new_element = " S " 
            new_atom_padded = f"{new_atom:<4}"
            new_line = line[:12] + new_atom_padded + line[16:-4] + new_element + "\n"
            new_lines[line_num - 1] = new_line

    with open('mol.pdb', 'w') as f:
        f.writelines(new_lines)
    
    with open('line_number.dat', 'w') as f:
        f.write('\n'.join(line_mappings))

def main():
    # Check if line_number.dat exists
    if not os.path.exists('line_number.dat'):
        # If it doesn't exist, process the PDB file
        process_pdb('COMPLEX.pdb')
    else:
        # If line_number.dat exists, process the MOL2 file
        atom_info = extract_atom_types_from_mol2('ONE.mol2')
        if not atom_info:
            print("Error processing MOL2 file")
            return
        # Then process COMPLEX.mol2 using both Reference_atom_type.dat and bond data from ONE.mol2
        if os.path.exists('Reference_atom_type.dat'):
            if replace_atom_types_and_bonds_from_mol2():
                print(" ")
            else:
                print("Error replacing atom types and bonds.")

if __name__ == "__main__":
    main()
