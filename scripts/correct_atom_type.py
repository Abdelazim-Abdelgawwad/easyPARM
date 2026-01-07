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
#                              |  $$$$$$/              Ver. 4.20 - 1 January 2026                                 #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


#Extract atom IDs and types from mol2 file.
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


#Read dihedral definitions from dihedral.dat file.
def read_dihedrals(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    dihedrals = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 4:
            try:
                atom1, atom2, atom3, atom4 = map(int, parts[:4])
                value = float(parts[4]) if len(parts) >= 5 else 0.0
                dihedrals.append((atom1, atom2, atom3, atom4, value))
            except ValueError:
                continue
    return dihedrals


#Normalize dihedral
def normalize_dihedral(atom1, atom2, atom3, atom4):
    forward = (atom1, atom2, atom3, atom4)
    reverse = (atom4, atom3, atom2, atom1)
    return min(forward, reverse)


#Check for cp cp cp cp or cq cq cq cq patterns and return modifications.
def check_and_modify_atom_types(atom_types, dihedrals):
    modifications = {}
    seen_dihedrals = set()
    
    for dihedral in dihedrals:
        atom1, atom2, atom3, atom4, value = dihedral
        
        normalized = normalize_dihedral(atom1, atom2, atom3, atom4)
        if normalized in seen_dihedrals:
            continue
        seen_dihedrals.add(normalized)
        
        type1 = atom_types.get(atom1)
        type2 = atom_types.get(atom2)
        type3 = atom_types.get(atom3)
        type4 = atom_types.get(atom4)
        
        if None in [type1, type2, type3, type4]:
            continue
        
        types = [type1, type2, type3, type4]
        
        if all(t == 'cp' for t in types):
            modifications[atom3] = 'cq'
            modifications[atom4] = 'cq'
        elif all(t == 'cq' for t in types):
            modifications[atom1] = 'cp'
            modifications[atom2] = 'cp'
    
    return modifications


#Verify no cp cp cp cp or cq cq cq cq patterns exist.
def verify_no_forbidden_patterns(atom_types, dihedrals):
    forbidden_patterns = []
    seen_dihedrals = set()
    
    for dihedral in dihedrals:
        atom1, atom2, atom3, atom4, value = dihedral
        
        normalized = normalize_dihedral(atom1, atom2, atom3, atom4)
        if normalized in seen_dihedrals:
            continue
        seen_dihedrals.add(normalized)
        
        types = [
            atom_types.get(atom1),
            atom_types.get(atom2),
            atom_types.get(atom3),
            atom_types.get(atom4)
        ]
        
        if None in types:
            continue
        
        if all(t == 'cp' for t in types):
            forbidden_patterns.append((dihedral, 'cp cp cp cp'))
        elif all(t == 'cq' for t in types):
            forbidden_patterns.append((dihedral, 'cq cq cq cq'))
    
    return forbidden_patterns


#Write modified mol2 file preserving exact format
def write_modified_mol2(input_file, output_file, modifications):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    is_atom_section = False
    modified_lines = []
    
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            is_atom_section = True
            modified_lines.append(line)
            continue
        
        if line.startswith("@<TRIPOS>BOND"):
            is_atom_section = False
            modified_lines.append(line)
            continue
        
        if is_atom_section:
            parts = line.split()
            if len(parts) >= 6:
                atom_id = int(parts[0])
                if atom_id in modifications:
                    old_type = parts[5]
                    new_type = modifications[atom_id]
                    
                    temp_parts = line.split()
                    if len(temp_parts) >= 6 and temp_parts[5] == old_type:
                        search_str = line
                        for i in range(5):
                            field_start = search_str.find(temp_parts[i])
                            search_str = search_str[field_start + len(temp_parts[i]):]
                        
                        type_pos = search_str.find(old_type)
                        actual_pos = len(line) - len(search_str) + type_pos
                        modified_line = line[:actual_pos] + new_type + line[actual_pos + len(old_type):]
                        modified_lines.append(modified_line)
                    else:
                        modified_lines.append(line)
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)
    
    with open(output_file, 'w') as file:
        file.writelines(modified_lines)


def main():
    mol2_file = "COMPLEX.mol2"
    dihedral_file = "dihedral.dat"
    output_file = "COMPLEX_modified_atom_type.mol2"
    
    # Read files
    atom_types = extract_atom_types_from_mol2(mol2_file)
    dihedrals = read_dihedrals(dihedral_file)
    # Check for patterns and get modifications
    modifications = check_and_modify_atom_types(atom_types, dihedrals)
    
    # Apply modifications
    for atom_id, new_type in modifications.items():
        atom_types[atom_id] = new_type
    # Verify no forbidden patterns remain
    forbidden = verify_no_forbidden_patterns(atom_types, dihedrals)
    # Write output
    write_modified_mol2(mol2_file, output_file, modifications)

if __name__ == "__main__":
    main()
