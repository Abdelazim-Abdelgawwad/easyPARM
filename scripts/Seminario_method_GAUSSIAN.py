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
#                              |  $$$$$$/              Ver. 3.00 - 12 January 2024                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import numpy as np
from collections import defaultdict
import sys

# function for gaussian input
def parse_gau(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None

    hessian_lines = []
    coord_lines = []
    charge_lines = []
    start_parsing_hessian = start_parsing_coords = start_parsing_charges = False
    coord_line_counter = 0
    
    for line in lines:
        # Parse Hessian
        if "Force constants in Cartesian coordinates" in line:
            start_parsing_hessian = True
            continue
        elif start_parsing_hessian and "Final forces over" in line:
            start_parsing_hessian = False
            continue
        elif start_parsing_hessian:
            hessian_lines.append(line.strip())
        
        # Parse coordinates
        if "Input orientation:" in line:
            start_parsing_coords = True
            coord_line_counter = 0
            coord_lines = []
            continue
        elif start_parsing_coords:
            coord_line_counter += 1
            if coord_line_counter > 4:
                if line.strip().startswith("-----"):
                    start_parsing_coords = False
                else:
                    coord_lines.append(line.strip())
        
        # Parse Mulliken charges
        if "Mulliken charges:" in line or "Mulliken atomic charges:" in line:
            start_parsing_charges = True
            charge_lines = []
            continue
        elif start_parsing_charges and ("Sum of Mulliken charges" in line or "Sum of Mulliken atomic charges" in line):
            start_parsing_charges = False
        elif start_parsing_charges:
            charge_lines.append(line.strip())
    
    coordinates, atomic_numbers = extract_coordinates(coord_lines)
    hessian_size = 3 * len(atomic_numbers)
    hessian_matrix = parse_blocked_hessian(hessian_lines, hessian_size)
    charges = extract_charges(charge_lines)
    
    return coordinates, hessian_matrix, atomic_numbers, charges

# function for checkpoint input
def parse_fchk(filename):
    def parse_array(lines, startline, endline):
        arr = []
        in_section = False
        for line in lines:
            if line.startswith(startline):
                in_section = True
                continue
            if in_section and line.startswith(endline):
                break
            if in_section:
                try:
                    arr.extend(map(float, line.split()))
                except ValueError:
                    break
        return np.array(arr, dtype=float)

    with open(filename, "r") as f:
        content = f.readlines()

    # Try different end markers for coordinates
    end_markers = ['Number of symbols in', 'Force Field', 'Atomic numbers']
    for end_marker in end_markers:
        crds = parse_array(content, 'Current cartesian coordinates', end_marker)
        if len(crds) > 0:
            break

    # Try different end markers for Hessian
    end_markers = ['Nonadiabatic coupling', 'Dipole Moment', 'Vibrational Anharmonic']
    for end_marker in end_markers:
        hess_lower = parse_array(content, 'Cartesian Force Constants', end_marker)
        if len(hess_lower) > 0:
            break

    atomic_numbers = parse_array(content, 'Atomic numbers', 'Nuclear charges')
    charges = parse_array(content, 'Mulliken Charges', 'Cartesian Gradient')
 
    if len(crds) == 0 or len(hess_lower) == 0 or len(atomic_numbers) == 0:
        raise ValueError("Failed to parse required data from the file.")

    natoms = len(crds) // 3
    dim = 3 * natoms
    hess = np.zeros((dim, dim))
    idx = 0
    for i in range(dim):
        for j in range(i+1):
            hess[i, j] = hess_lower[idx]
            hess[j, i] = hess_lower[idx]
            idx += 1
    
    return crds.reshape(-1, 3), hess, atomic_numbers.astype(int), charges

# Helper functions for parse_gau
def extract_coordinates(coord_lines):
    coordinates = []
    atomic_numbers = []
    
    for line in coord_lines:
        parts = line.split()
        if len(parts) >= 6:
            atomic_numbers.append(int(parts[1]))
            coordinates.append(list(map(float, parts[3:6])))
    
    return np.array(coordinates), np.array(atomic_numbers)

def extract_charges(charge_lines):
    charges = []
    for line in charge_lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                charges.append(float(parts[-1]))
            except ValueError:
                continue
    return np.array(charges)

def parse_blocked_hessian(hessian_lines, hessian_size):
    hessian_matrix = np.zeros((hessian_size, hessian_size))
    current_col_start = 0
    
    for line in hessian_lines:
        parts = line.split()
        
        if all(p.isdigit() for p in parts):
            current_col_start = int(parts[0]) - 1
            continue
        
        row = int(parts[0]) - 1
        values = [float(val.replace('D', 'E')) for val in parts[1:]]
        
        for i, value in enumerate(values):
            col = current_col_start + i
            if col < hessian_size and row < hessian_size:
                hessian_matrix[row, col] = value
                hessian_matrix[col, row] = value
    
    return hessian_matrix

BOHR_TO_ANGSTROM = 0.529177
HARTREE_TO_KCAL_MOL = 627.509474

#Calculate distance between two coordinates.
def calculate_distance(coord1, coord2, is_fchk=False):

    distance = np.linalg.norm(coord1 - coord2)
    if is_fchk:
        return distance * BOHR_TO_ANGSTROM
    return distance

def extract_sub_hessian(hessian, i, j):
    indices = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
    sub_hess = - hessian[np.ix_(indices[:3], indices[3:])]
    return sub_hess

def calculate_bond_force_constant(hessian, coordinates, idx1, idx2):
    # extract sub-hessians for bond involved in the angle
    sub_hessian = extract_sub_hessian(hessian, idx1, idx2)

    # Calculate the vector from atom1 to atom2
    vec12 = coordinates[idx2] - coordinates[idx1]
    vec12 = vec12 
    # Normalize the vector to get the unit vector along the bond
    unit_vec12 = vec12 / np.linalg.norm(vec12)

    # Calculate eigenvalues and eigenvectors of the sub-Hessian
    eigvals, eigvecs = np.linalg.eigh(sub_hessian)

    # Calculate the force constant using the Seminario Method
    # This involves projecting each eigenvector onto the bond unit vector
    # and weighting it by the corresponding eigenvalue
    force_constant = sum((np.dot(eigvecs[:, i], unit_vec12)**2) * eigvals[i] for i in range(3))

    # Convert force constant from atomic units (Hartree/Bohr^2) to kcal/mol/A^2
    force_constant_in_kcal = 627.509474 * force_constant # 627.509474 is the conversion factor

    # Apply the harmonic approximation (factor of 2)
    force_constant_harmonic = 2 * force_constant_in_kcal

    return force_constant_harmonic

def seminario_method(hessian, coordinates, atom_pairs, atomic_numbers):
    bonds = []
    for atom1, atom2 in atom_pairs:
        idx1, idx2 = atom1 - 1, atom2 - 1
        force_constant = calculate_bond_force_constant(hessian, coordinates, idx1, idx2)
        bonds.append((atom1, atom2, force_constant))
    return bonds

def read_distances(filename):

    atom_pairs = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                parts = line.split()
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                atom_pairs.append((atom1, atom2))
    return atom_pairs

def read_angles(filename):
    angles = []
    with open(filename, "r") as f:
        for line in f:
            if line.strip():
                parts = line.split()
                atom1 = int(parts[0])
                atom2 = int(parts[1])
                atom3 = int(parts[2])
                angle_value = float(parts[3])
                angles.append((atom1, atom2, atom3, angle_value))
    return angles

def calculate_angle(coord1, coord2, coord3):
    v1 = coord1 - coord2
    v2 = coord3 - coord2

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))

    return np.degrees(angle)

#Calculate angle force constant with proper unit handling for both Gaussian output and fchk files.    
def calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3, is_fchk=False):
    # Extract sub-hessians for both bonds involved in the angle
    sub_hessian_AB = extract_sub_hessian(hessian, idx1, idx2)
    sub_hessian_CB = extract_sub_hessian(hessian, idx3, idx2)

    # Calculate vectors
    vec_AB = coordinates[idx1] - coordinates[idx2]
    vec_CB = coordinates[idx3] - coordinates[idx2]
    
    # Convert vectors to Bohr if they're in Angstroms (Gaussian output)
    BOHR_TO_ANGSTROM = 0.529177
    HARTREE_TO_KCAL_MOL = 627.509474
    
    if not is_fchk:
        vec_AB = vec_AB / BOHR_TO_ANGSTROM 
        vec_CB = vec_CB / BOHR_TO_ANGSTROM

    # Calculate unit vectors
    u_AB = vec_AB / np.linalg.norm(vec_AB)
    u_CB = vec_CB / np.linalg.norm(vec_CB)

    # Calculate u_N (equation 11)
    u_N = np.cross(u_CB, u_AB) / np.linalg.norm(np.cross(u_CB, u_AB))

    # Calculate u_PA and u_PC (equations 12 and 13)
    u_PA = np.cross(u_N, u_AB)
    u_PC = np.cross(u_CB, u_N)

    # Calculate eigenvalues and eigenvectors
    eigvals_AB, eigvecs_AB = np.linalg.eigh(sub_hessian_AB)
    eigvals_CB, eigvecs_CB = np.linalg.eigh(sub_hessian_CB)

    # Calculate R_AB and R_CB (now in Bohr)
    R_AB = np.linalg.norm(vec_AB)
    R_CB = np.linalg.norm(vec_CB)

    # Implement equation (14)
    k_theta_AB = sum((np.dot(eigvecs_AB[:, i], u_PA)**2) * eigvals_AB[i] for i in range(3))
    k_theta_CB = sum((np.dot(eigvecs_CB[:, i], u_PC)**2) * eigvals_CB[i] for i in range(3))

    k_theta = 1 / (1 / (R_AB**2 * k_theta_AB) + 1 / (R_CB**2 * k_theta_CB))
    k_theta = abs(k_theta)

    # Convert from Hartree/Bohr²/radian² to kcal/mol/radian²
    k_theta_kcal_rad = 2 * k_theta * HARTREE_TO_KCAL_MOL * ( BOHR_TO_ANGSTROM ** 2)
    
    # For fchk files, we need the additional Bohr to Angstrom conversion

    return k_theta_kcal_rad

def read_dihedrals(filename):
    dihedrals = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                if len(parts) == 5:
                    atom1, atom2, atom3, atom4, dihedral_value = map(float, parts)
                    dihedrals.append((int(atom1), int(atom2), int(atom3), int(atom4), dihedral_value))
    return dihedrals

# Determine connectivity between atoms based on distance
def get_connectivity(coordinates, atomic_numbers, bond_threshold=1.7):
    connectivity = defaultdict(list)
    n_atoms = len(atomic_numbers)
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            if distance <= bond_threshold:
                connectivity[i].append(j)
                connectivity[j].append(i)
    return connectivity

# Find ring structures in the molecule.
def find_rings(connectivity):
    rings = []
    def dfs(node, path):
        if len(path) > 2 and node == path[0]:
            rings.append(tuple(sorted(path)))
            return
        for neighbor in connectivity[node]:
            if neighbor not in path[1:]:
                dfs(neighbor, path + [neighbor])
    
    for start_node in connectivity:
        dfs(start_node, [start_node])
    
    return list(set(rings))

# Generate a unique signature for an atom based on its properties and environment.
def get_atom_signature(atom, atomic_numbers, connectivity, rings, charges, charge_tolerance=0.01):
    atom_type = atomic_numbers[atom]
    charge = charges[atom]
    neighbors = tuple(sorted([atomic_numbers[n] for n in connectivity[atom]]))
    ring_types = tuple(sorted(tuple(sorted(atomic_numbers[i] for i in ring)) for ring in rings if atom in ring))
    return (atom_type, round(charge / charge_tolerance) * charge_tolerance, neighbors, ring_types)

# Detect atoms with similar environments and charges, considering only the same atom type.
def detect_similar_atoms(coordinates, atomic_numbers, charges):
    connectivity = get_connectivity(coordinates, atomic_numbers)
    rings = find_rings(connectivity)
    signatures = {}
    similar_atoms = defaultdict(list)

    for i in range(len(atomic_numbers)):
        sig = get_atom_signature(i, atomic_numbers, connectivity, rings, charges)
        atom_type = atomic_numbers[i]
        if sig in signatures and atomic_numbers[signatures[sig]] == atom_type:
            similar_atoms[signatures[sig]].append(i)
        else:
            signatures[sig] = i
            similar_atoms[i] = [i]

    return similar_atoms

# Write similar atoms to a file, considering only the same atom type.
def write_similar_atoms(similar_atoms, atomic_numbers, charges, filename="similar.dat"):
    with open(filename, "w") as similar_file:
        for reference, group in similar_atoms.items():
            if len(group) > 1:
                for atom in group:
                    if atom != reference and atomic_numbers[atom] == atomic_numbers[reference]:
                        similar_file.write(f"{atom+1:5d} {reference+1:5d}\n")

def main(input_file, file_type):
    distance_file = "distance.dat"
    angle_file = "angle.dat"
    dihedral_file = "dihedral.dat"

    # Choose parsing method based on file_type
    is_fchk = file_type in [3, 4]
    if file_type == 2:
        coordinates, hessian, atomic_numbers, charges = parse_gau(input_file)
    else:
        coordinates, hessian, atomic_numbers, charges = parse_fchk(input_file)

    connectivity = get_connectivity(coordinates, atomic_numbers)
    atom_pairs = read_distances(distance_file)
    angle_definitions = read_angles(angle_file)
    dihedral_definitions = read_dihedrals(dihedral_file)

    bonds = seminario_method(hessian, coordinates, atom_pairs, atomic_numbers)
    similar_atoms = detect_similar_atoms(coordinates, atomic_numbers, charges)
    write_similar_atoms(similar_atoms, atomic_numbers, charges)

    with open("bond_angle_dihedral_data.dat", "w") as file:
        # Write bond information
        file.write("\nBond Information:\n")
        for bond in bonds:
            atom1, atom2 = bond[0] - 1, bond[1] - 1
            calculated_distance = calculate_distance(coordinates[atom1], coordinates[atom2], is_fchk)
            force_constant = bond[2]
            file.write(f"{bond[0]:5d} {bond[1]:5d} {calculated_distance:10.3f} {force_constant:15.2f}\n")

        # Write angle information
        file.write("\nAngle Information:\n")
        for angle_def in angle_definitions:
            atom1, atom2, atom3, _ = angle_def
            idx1, idx2, idx3 = atom1 - 1, atom2 - 1, atom3 - 1
            calculated_angle = calculate_angle(coordinates[idx1], coordinates[idx2], coordinates[idx3])
            force_constant = calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3, is_fchk)
            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {calculated_angle:10.3f} {force_constant:15.2f}\n")

        # Write dihedral information
        file.write("\nDihedral Information:\n")
        for dihedral_def in dihedral_definitions:
            atom1, atom2, atom3, atom4, dihedral_value = dihedral_def
            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {atom4:5d} {dihedral_value:10.3f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py input_file file_type")
        print("file_type: 2 for Gaussian output file, 3 or 4 for Gaussian fchk file")
        sys.exit(1)

    input_file = sys.argv[1]
    file_type = int(sys.argv[2])

    if file_type not in [2, 3, 4]:
        print("Error: file_type must be 2 (Gaussian output) or 3/4 (Gaussian fchk)")
        sys.exit(1)

    main(input_file, file_type)
