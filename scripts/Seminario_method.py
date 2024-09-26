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
#                              |  $$$$$$/              Ver. 1.00 - 26 September 2024                              #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, ValÃ¨ncia 46071, Spain          #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de ValÃ¨ncia. E-mail: abdelazim.abdelgawwad@uv.es        #
###################################################################################################################


import numpy as np
from collections import defaultdict
import sys

# Read the Gaussian fchk file to extract coordinates, Hessian, atomic numbers, and charges.
def parse_fchk(filename):
    def parse_array(f, startline, endline):
        arr = []
        in_section = False
        for line in f:
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
        content = f.read()

    # Try different end markers for coordinates
    end_markers = ['Number of symbols in', 'Force Field', 'Atomic numbers']
    for end_marker in end_markers:
        crds = parse_array(content.splitlines(), 'Current cartesian coordinates', end_marker)
        if len(crds) > 0:
            break

    # Try different end markers for Hessian
    end_markers = ['Nonadiabatic coupling', 'Dipole Moment', 'Vibrational Anharmonic']
    for end_marker in end_markers:
        hess_lower = parse_array(content.splitlines(), 'Cartesian Force Constants', end_marker)
        if len(hess_lower) > 0:
            break
    # Extract the atomic number and charges
    atomic_numbers = parse_array(content.splitlines(), 'Atomic numbers', 'Nuclear charges')
    charges = parse_array(content.splitlines(), 'Mulliken Charges', 'Cartesian Gradient')
 
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

# Extract the sub hessian matrix (3,3)
def extract_sub_hessian(hessian, i, j):
    indices = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2]
    sub_hess = - hessian[np.ix_(indices[:3], indices[3:])]
    return sub_hess


# Extract the 3x3 sub-Hessian for the bond between atoms idx1 and idx2
def calculate_bond_force_constant(hessian, coordinates, idx1, idx2):
    sub_hessian = extract_sub_hessian(hessian, idx1, idx2)

    # Calculate the vector from atom1 to atom2
    vec12 = coordinates[idx2] - coordinates[idx1]

    # Normalize the vector to get the unit vector along the bond
    unit_vec12 = vec12 / np.linalg.norm(vec12)

    # Calculate eigenvalues and eigenvectors of the sub-Hessian
    eigvals, eigvecs = np.linalg.eigh(sub_hessian)

    # Calculate the force constant using the Seminario Method
    # This involves projecting each eigenvector onto the bond unit vector
    # and weighting it by the corresponding eigenvalue
    force_constant = sum((np.dot(eigvecs[:, i], unit_vec12)**2) * eigvals[i] for i in range(3))

    # Convert force constant from atomic units (Hartree/Bohr^2) to kcal/mol/Ãƒâ€¦^2
    force_constant_in_kcal = 627.509474 * force_constant  # 627.509474 is the conversion factor

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

# Read the distance that contain the pairwise atoms info
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

def calculate_distance(coord1, coord2):
    bohr_to_angstrom = 0.529177
    return np.linalg.norm(coord1 - coord2) * bohr_to_angstrom

# Read the angle file that the three atoms that form the angle
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
    
#Calculate the force constant for an angle between atom1-atom2-atom3 using the Seminario method.    
#Return the force constant in kcal/mol/radian².
def calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3):
    # Extract sub-hessians
    sub_hessian_AB = extract_sub_hessian(hessian, idx1, idx2)
    sub_hessian_CB = extract_sub_hessian(hessian, idx3, idx2)

    # Calculate unit vectors
    vec_AB = coordinates[idx1] - coordinates[idx2]
    vec_CB = coordinates[idx3] - coordinates[idx2]
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

    # Calculate R_AB and R_CB (in Bohr)
    R_AB = np.linalg.norm(vec_AB)
    R_CB = np.linalg.norm(vec_CB)

    # Implement equation (14)
    k_theta_AB = sum((np.dot(eigvecs_AB[:, i], u_PA)**2) * eigvals_AB[i] for i in range(3))
    k_theta_CB = sum((np.dot(eigvecs_CB[:, i], u_PC)**2) * eigvals_CB[i] for i in range(3))

    k_theta = 1 / (1 / (R_AB**2 * k_theta_AB) + 1 / (R_CB**2 * k_theta_CB))
    k_theta = abs(k_theta)

    # Convert from Hartree/Bohr²/radian² to kcal/mol/radian²
    hartree_to_kcal_mol = 627.509474
    bohr_to_angstrom = 0.529177
    k_theta_kcal_rad = 2 * k_theta * hartree_to_kcal_mol * (bohr_to_angstrom ** 2)

    return k_theta_kcal_rad

# Read the dihedral info from file
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

#Determine connectivity between atoms based on distance.
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

#Find ring structures in the molecule.
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

#Generate a unique signature for an atom based on its properties and environment
def get_atom_signature(atom, atomic_numbers, connectivity, rings, charges, charge_tolerance=0.015):
    atom_type = atomic_numbers[atom]
    charge = charges[atom]
    neighbors = tuple(sorted([atomic_numbers[n] for n in connectivity[atom]]))
    ring_types = tuple(sorted(tuple(sorted(atomic_numbers[i] for i in ring)) for ring in rings if atom in ring))
    return (atom_type, round(charge / charge_tolerance) * charge_tolerance, neighbors, ring_types)

#Detect atoms with similar environments in the molecule.
def detect_similar_atoms(coordinates, atomic_numbers, charges):
    connectivity = get_connectivity(coordinates, atomic_numbers)
    rings = find_rings(connectivity)
    signatures = {}
    similar_atoms = defaultdict(list)

    for i in range(len(atomic_numbers)):
        sig = get_atom_signature(i, atomic_numbers, connectivity, rings, charges)
        if sig in signatures:
            similar_atoms[signatures[sig]].append(i)
        else:
            signatures[sig] = i
            similar_atoms[i] = [i]

    return similar_atoms

#Write similar atoms to a file.
def write_similar_atoms(similar_atoms, atomic_numbers, charges, filename="similar.dat"):
    with open(filename, "w") as similar_file:
        for reference, group in similar_atoms.items():
            if len(group) > 1:
                for atom in group:
                    if atom != reference:
                        similar_file.write(f"{atom+1:5d} {reference+1:5d}\n")

# Run all the function
def main(fchk_file):
    # File paths for storing computed distance, angle, and dihedral data
    distance_file = "distance.dat"
    angle_file = "angle.dat"
    dihedral_file = "dihedral.dat"
    
    # Read data from the formatted checkpoint file (fchk), extracting coordinates, Hessian, atomic numbers, and charges
    coordinates, hessian, atomic_numbers, charges = parse_fchk(fchk_file)
    connectivity = get_connectivity(coordinates, atomic_numbers)

    atom_pairs = read_distances(distance_file)
    angle_definitions = read_angles(angle_file)
    dihedral_definitions = read_dihedrals(dihedral_file)
    
    # Apply the Seminario method to calculate bond force constants using Hessian, coordinates, and atomic pairs
    bonds = seminario_method(hessian, coordinates, atom_pairs, atomic_numbers)

    # Detect atoms with similar properties based on their coordinates, atomic numbers, and charges
    similar_atoms = detect_similar_atoms(coordinates, atomic_numbers, charges)
    write_similar_atoms(similar_atoms, atomic_numbers, charges)

    # Open an output file to store bond, angle, and dihedral information
    with open("bond_angle_dihedral_data.dat", "w") as file:
        file.write("\nBond Information:\n")
        for bond in bonds:
            atom1, atom2 = bond[0] - 1, bond[1] - 1
            calculated_distance = calculate_distance(coordinates[atom1], coordinates[atom2])
            force_constant = bond[2]
            file.write(f"{bond[0]:5d} {bond[1]:5d} {calculated_distance:10.3f} {force_constant:15.2f}\n")

        file.write("\nAngle Information:\n")
        for angle_def in angle_definitions:
            atom1, atom2, atom3, _ = angle_def
            idx1, idx2, idx3 = atom1 - 1, atom2 - 1, atom3 - 1  # Convert to 0-based indexing

            # Calculate the angle
            calculated_angle = calculate_angle(coordinates[idx1], coordinates[idx2], coordinates[idx3])

            force_constant_rad = calculate_angle_force_constant(hessian, coordinates, idx1, idx2, idx3)

            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {calculated_angle:10.3f} {force_constant_rad:15.2f}\n")

        file.write("\nDihedral Information:\n")
        for dihedral_def in dihedral_definitions:
            atom1, atom2, atom3, atom4, dihedral_value = dihedral_def
            file.write(f"{atom1:5d} {atom2:5d} {atom3:5d} {atom4:5d} {dihedral_value:10.3f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py path/to/file.fchk")
        sys.exit(1)

    main(sys.argv[1])

