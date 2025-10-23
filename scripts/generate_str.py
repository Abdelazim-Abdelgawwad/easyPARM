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

import sys
import re

# Previous atomic_masses dictionary and helper functions remain the same...
atomic_masses = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.63,
    "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 97.91, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33,
    "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49,
    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237, "Pu": 244,
    "Am": 243, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252, "Fm": 257, "Md": 258, "No": 259, "Lr": 266,
    "Rf": 267, "Db": 270, "Sg": 271, "Bh": 270, "Hs": 277, "Mt": 276, "Ds": 281, "Rg": 282, "Cn": 285,
    "Nh": 286, "Fl": 289, "Mc": 290, "Lv": 293, "Ts": 294, "Og": 294
}

def read_forcefield_data(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        bonds = []
        angles = []
        reading_bonds = False
        reading_angles = False
        
        for line in lines:
            line = line.strip()
            
            if line == "BOND":
                reading_bonds = True
                reading_angles = False
                continue
            elif line == "ANGLE":
                reading_bonds = False
                reading_angles = True
                continue
            elif not line:  # Skip empty lines
                continue
            
            if reading_bonds and "-" in line:
                parts = line.split()
                if parts:
                    atom_pair = parts[0].split("-")
                    if len(atom_pair) == 2:
                        bonds.append((atom_pair[0], atom_pair[1]))
            
            elif reading_angles and "-" in line:
                parts = line.split()
                if parts:
                    atoms = parts[0].split("-")
                    if len(atoms) == 3:
                        angles.append((atoms[0], atoms[1], atoms[2]))
        
        return bonds, angles
    
    except Exception as e:
        print(f"Error reading forcefield file {file_path}: {str(e)}")
        sys.exit(1)

def extract_atoms_from_mol2(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        atom_data = []
        total_charge = 0.0
        is_atom_section = False

        for line in lines:
            line = line.strip()
            
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            elif line.startswith("@<TRIPOS>"):
                is_atom_section = False
                continue

            if is_atom_section and line:
                parts = line.split()
                if len(parts) >= 9:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    atom_type = parts[5]
                    resid_id = parts[7]
                    charge = float(parts[8])
                    
                    element_symbol = re.sub(r'[0-9]', '', atom_name)  # Remove numbers
                    # Capitalize first letter for standardization
                    element_symbol = element_symbol.capitalize()
                    
                    
                    #if atom_type == "RU":
                    #    element_symbol = "Ru"
                    #else:
                    #    element_symbol = "C" if atom_type.startswith("C") else "H"
                    
                    atomic_mass = atomic_masses.get(element_symbol, 0.0)
                    
                    total_charge += charge
                    atom_data.append((atom_id, atom_name, atom_type, atomic_mass, element_symbol, charge))

        return atom_data, total_charge, resid_id

    except Exception as e:
        print(f"Error reading mol2 file {file_path}: {str(e)}")
        sys.exit(1)

def print_complete_output(atom_data, total_charge, bonds, angles, resid_id): 
    if not atom_data:
        print("ERROR: No atom data found in the file!")
        return

    # MASS Section
    print("* >>>>>>>>>>>>>>>>> MOL TOPOLOGY USING easyPARM")
    print("* (FEB 2025)")
    print("\n")
    
    unique_atom_types = set()
    for atom_id, _, atom_type, atomic_mass, element_symbol, _ in atom_data:
        if atom_type not in unique_atom_types:
            print(f"MASS   {atom_id:<3} {atom_type:<5} {atomic_mass:.4f}   {element_symbol:<2}  ")
            unique_atom_types.add(atom_type)

    # RESI Section
    print(f"\nRESI {resid_id:<1}       {total_charge:.3f}")
    print("GROUP")

    # ATOM Section
    for _, atom_name, atom_type, _, _, charge in atom_data:
        print(f"ATOM   {atom_name:<4} {atom_type:<5} {charge:>9.6f}")

    # Blank line then BOND section
    print("\n")
    for atom1, atom2 in bonds:
        print(f"BOND  {atom1}  {atom2}")

    # Blank line then ANGLE section
    print("\n")
    for atom1, atom2, atom3 in angles:
        print(f"ANGLE  {atom1}  {atom2}  {atom3}")
    print("\n")
    print("END")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_str.py <mol2_file> <forcefield_file>")
        sys.exit(1)

    mol2_file = sys.argv[1]
    forcefield_file = sys.argv[2]
    
    # Read both files
    atom_data, total_charge, resid_id = extract_atoms_from_mol2(mol2_file)
    bonds, angles = read_forcefield_data(forcefield_file)
    
    # Print complete output
    print_complete_output(atom_data, total_charge, bonds, angles, resid_id)
