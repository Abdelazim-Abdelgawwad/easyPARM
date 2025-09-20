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
#                              |  $$$$$$/              Ver. 4.10 - 20 September 2025                              #
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
import periodictable

def extract_coordinates_from_file(file_path):
    coordinates = []
    atom_names = []
    reading_atoms = False
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == '$atoms':
                reading_atoms = True
                num_atoms = int(next(file).strip())
                continue
                
            if reading_atoms:
                if line.strip().startswith('$'):
                    break
                if not line.strip():
                    continue
                    
                parts = line.split()
                if len(parts) >= 5:
                    atom_names.append(parts[0].strip())
                    coordinates.append([float(x) for x in parts[-3:]])
    
    coordinates = np.array(coordinates)
    atomic_numbers = np.array([periodictable.elements.symbol(name).number for name in atom_names])
    
    return coordinates, atomic_numbers, atom_names

def extract_hessian_from_file(file_path):
    hessian_data = []
    reading_hessian = False
    hessian_size = 0
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    for i, line in enumerate(lines):
        if line.strip().startswith('$hessian'):
            reading_hessian = True
            hessian_size = int(lines[i+1].strip())
            start_index = i + 2
            break
    
    if not reading_hessian:
        raise ValueError("Hessian section not found in the file.")
    
    n_atoms = hessian_size // 3
    expected_elements = hessian_size * hessian_size  # Full matrix
    
    hessian_matrix = np.zeros((hessian_size, hessian_size))
    current_col = 0
    
    while current_col < hessian_size:
        try:
            # Skip the line with column indices
            start_index += 1
            if start_index >= len(lines):
                print(f"Reached end of file unexpectedly. Current column: {current_col}")
                break
            
            num_columns = min(5, hessian_size - current_col)
            end_col = current_col + num_columns
            
            for row in range(hessian_size):
                line_index = start_index + row
                if line_index >= len(lines):
                    break
                line = lines[line_index].strip()
                if not line:  # Skip empty lines
                    continue
                parts = line.split()
                if len(parts) <= 1:  # Skip lines with only an index
                    continue
                parts = parts[1:]  # Skip the row index
                for col, value in enumerate(parts):
                    global_col = current_col + col
                    if global_col >= end_col:
                        break
                    try:
                        hessian_matrix[row, global_col] = float(value.replace('D', 'E').replace('E', 'e'))
                    except ValueError:
                        print(f"Error converting value to float: '{value}' at row {row}, col {global_col}")
            
            current_col = end_col
            start_index += hessian_size  # Move to the next block
            # Skip any blank lines
            while start_index < len(lines) and not lines[start_index].strip():
                start_index += 1
        except Exception as e:
            break

    return hessian_matrix

def extract_mulliken_charges(file_path):
    charges = []
    reading_charges = False
    expected_atom_count = None
    
    with open(file_path, 'r') as file:
        for line in file:
            if "MULLIKEN ATOMIC CHARGES" in line:
                reading_charges = True
                # Try to extract the number of atoms from the next line
                next_line = next(file, '').strip()
                if next_line.startswith('Number of atoms'):
                    expected_atom_count = int(next_line.split()[-1])
                continue
            
            if reading_charges:
                if line.strip() == '':
                    continue
                if "Sum of atomic charges" in line:
                    break
                
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        index = int(parts[0])
                        atom = parts[1]
                        charge = float(parts[-1])  # Assume charge is the last column
                        charges.append(charge)
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing line: {line.strip()}")
                        print(f"Error details: {e}")
                        continue
    

    return np.array(charges)

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
    force_constant = abs(force_constant)
    # Convert force constant from atomic units (Hartree/Bohr^2) to kcal/mol/Ã^2
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

# Enhanced molecular connectivity analysis
def analyze_molecular_connectivity(coordinates, atomic_numbers, atom_names):
    
    bohr_to_angstrom = 0.529177
    natom = len(atomic_numbers)
    
    # Calculate distance matrix
    distance_matrix = np.zeros((natom, natom))
    for i in range(natom):
        for j in range(natom):
            if i != j:
                distance_matrix[i][j] = np.linalg.norm(coordinates[i] - coordinates[j]) * bohr_to_angstrom
    
    # Define typical bond lengths for connectivity (in Angstroms)
    bond_cutoffs = {
        ('H', 'C'): 1.2, ('H', 'N'): 1.2, ('H', 'O'): 1.2,
        ('C', 'C'): 1.8, ('C', 'N'): 1.8, ('C', 'O'): 1.8,
        ('N', 'N'): 1.8, ('N', 'O'): 1.8, ('O', 'O'): 1.8,
        ('H', 'S'): 1.5, ('C', 'S'): 2.0, ('S', 'S'): 2.3,
        ('H', 'P'): 1.5, ('C', 'P'): 2.0, ('P', 'P'): 2.3,
        # Metal bonds (more generous cutoffs)
        ('Metal', 'H'): 2.0, ('Metal', 'C'): 2.5, ('Metal', 'N'): 2.5,
        ('Metal', 'O'): 2.5, ('Metal', 'S'): 2.8, ('Metal', 'P'): 2.8,
        ('Metal', 'Metal'): 3.0
    }
    
    # Metal atoms set
    metal_atoms = {
        'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
        'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
        'LA', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
        'AC', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN'
    }
    
    def get_cutoff(elem1, elem2):
        elem1_is_metal = periodictable.elements[elem1].symbol.upper() in metal_atoms
        elem2_is_metal = periodictable.elements[elem2].symbol.upper() in metal_atoms
        
        if elem1_is_metal or elem2_is_metal:
            if elem1_is_metal and elem2_is_metal:
                return bond_cutoffs.get(('Metal', 'Metal'), 3.0)
            else:
                non_metal = elem2 if elem1_is_metal else elem1
                non_metal_symbol = periodictable.elements[non_metal].symbol.upper()
                return bond_cutoffs.get(('Metal', non_metal_symbol), 
                                      bond_cutoffs.get((non_metal_symbol, 'Metal'), 2.5))
        else:
            elem1_symbol = periodictable.elements[elem1].symbol.upper()
            elem2_symbol = periodictable.elements[elem2].symbol.upper()
            pair1 = (elem1_symbol, elem2_symbol)
            pair2 = (elem2_symbol, elem1_symbol)
            return bond_cutoffs.get(pair1, bond_cutoffs.get(pair2, 1.8))
    
    # Build connectivity matrix
    connectivity = np.zeros((natom, natom), dtype=bool)
    for i in range(natom):
        for j in range(i+1, natom):
            elem_i = atomic_numbers[i]
            elem_j = atomic_numbers[j]
            cutoff = get_cutoff(elem_i, elem_j)
            
            if distance_matrix[i][j] <= cutoff:
                connectivity[i][j] = True
                connectivity[j][i] = True
    
    # Create neighbor lists
    neighbors = {}
    for i in range(natom):
        neighbors[i] = []
        for j in range(natom):
            if connectivity[i][j]:
                neighbors[i].append(j)
    
    return connectivity, neighbors, distance_matrix

#Create a detailed chemical fingerprint for an atom based on its environment
def get_chemical_environment_fingerprint(atomic_numbers, atom_idx, neighbors, max_depth=3):
    
    def explore_environment(current_atom, depth, visited, path):
        if depth > max_depth:
            return []
        
        visited.add(current_atom)
        environment_info = []
        
        # Get immediate neighbors
        current_neighbors = []
        for neighbor_idx in neighbors[current_atom]:
            if neighbor_idx not in visited:
                neighbor_element = periodictable.elements[atomic_numbers[neighbor_idx]].symbol.upper()
                current_neighbors.append((neighbor_element, neighbor_idx))
        
        # Sort neighbors for consistent fingerprinting
        current_neighbors.sort(key=lambda x: (x[0], x[1]))
        
        for neighbor_element, neighbor_idx in current_neighbors:
            # Create path signature
            new_path = path + [neighbor_element]
            environment_info.append(tuple(new_path))
            
            # Recursively explore deeper
            if depth < max_depth:
                deeper_env = explore_environment(neighbor_idx, depth + 1, visited.copy(), new_path)
                environment_info.extend(deeper_env)
        
        return environment_info
    
    # Start exploration from the atom
    fingerprint = explore_environment(atom_idx, 0, set(), [periodictable.elements[atomic_numbers[atom_idx]].symbol.upper()])
    return tuple(sorted(fingerprint))

#Identify aromatic rings in the molecule using simple ring detection
def identify_aromatic_rings(atomic_numbers, connectivity, neighbors):
    
    natom = len(atomic_numbers)
    aromatic_atoms = set()
    
    # Find rings using DFS
    def find_rings_from_atom(start_atom, max_ring_size=8):
        def dfs(current, path, visited_edges):
            if len(path) > max_ring_size:
                return []
            
            rings_found = []
            for neighbor in neighbors[current]:
                edge = tuple(sorted([current, neighbor]))
                
                if neighbor == start_atom and len(path) >= 3:
                    # Found a ring
                    rings_found.append(path[:])
                elif neighbor not in path and edge not in visited_edges:
                    new_visited = visited_edges.copy()
                    new_visited.add(edge)
                    rings_found.extend(dfs(neighbor, path + [neighbor], new_visited))
            
            return rings_found
        
        return dfs(start_atom, [start_atom], set())
    
    # Find all rings
    all_rings = []
    processed = set()
    
    for atom in range(natom):
        if atom not in processed:
            rings = find_rings_from_atom(atom)
            for ring in rings:
                if len(ring) >= 5 and len(ring) <= 8:  # Reasonable ring sizes
                    ring_tuple = tuple(sorted(ring))
                    if ring_tuple not in [tuple(sorted(r)) for r in all_rings]:
                        all_rings.append(ring)
                        processed.update(ring)
    
    # Check for aromaticity (simplified heuristic)
    for ring in all_rings:
        if len(ring) in [5, 6]:  # Common aromatic ring sizes
            # Check if ring is mostly carbons with some sp2 character
            carbon_count = sum(1 for atom in ring if periodictable.elements[atomic_numbers[atom]].symbol.upper() == 'C')
            nitrogen_count = sum(1 for atom in ring if periodictable.elements[atomic_numbers[atom]].symbol.upper() == 'N')
            
            # Heuristic: if mostly C and N, likely aromatic
            if (carbon_count + nitrogen_count) >= len(ring) * 0.8:
                aromatic_atoms.update(ring)
    
    return aromatic_atoms

#Estimate carbon hybridization based on bonding pattern.
def analyze_carbon_hybridization(atomic_numbers, carbon_idx, neighbors):
    carbon_neighbors = neighbors[carbon_idx]
    num_neighbors = len(carbon_neighbors)
    
    # Count different atom types
    neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in carbon_neighbors]
    hydrogen_count = neighbor_elements.count('H')
    heavy_atom_count = num_neighbors - hydrogen_count
    
    if num_neighbors == 4:
        return 'sp3'
    elif num_neighbors == 3:
        return 'sp2'
    elif num_neighbors == 2:
        return 'sp'
    else:
        return 'unknown'

#Enhanced hydrogen environment analysis
def analyze_enhanced_hydrogen_environments(atomic_numbers, neighbors, connectivity):
    natom = len(atomic_numbers)
    
    # First identify aromatic rings
    aromatic_atoms = identify_aromatic_rings(atomic_numbers, connectivity, neighbors)
    
    # Categorize hydrogens with detailed environment analysis
    hydrogen_environments = {}
    
    for i in range(natom):
        if periodictable.elements[atomic_numbers[i]].symbol.upper() != 'H':
            continue
        
        # Get the heavy atom this H is bonded to
        heavy_neighbors = []
        for neighbor_idx in neighbors[i]:
            neighbor_element = periodictable.elements[atomic_numbers[neighbor_idx]].symbol.upper()
            if neighbor_element != 'H':
                heavy_neighbors.append((neighbor_idx, neighbor_element))
        
        if len(heavy_neighbors) != 1:
            continue  # Skip unusual cases
            
        heavy_atom_idx, heavy_element = heavy_neighbors[0]
        
        # Create detailed environment fingerprint
        env_fingerprint = get_chemical_environment_fingerprint(atomic_numbers, i, neighbors, max_depth=3)
        
        # Analyze specific environments
        if heavy_element == 'C':
            carbon_hybridization = analyze_carbon_hybridization(atomic_numbers, heavy_atom_idx, neighbors)
            is_aromatic_carbon = heavy_atom_idx in aromatic_atoms
            
            # Get carbon's neighbors (excluding this hydrogen)
            carbon_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            carbon_neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in carbon_neighbors]
            
            # Count different types of neighbors
            h_on_carbon = carbon_neighbor_elements.count('H')
            c_on_carbon = carbon_neighbor_elements.count('C')
            n_on_carbon = carbon_neighbor_elements.count('N')
            o_on_carbon = carbon_neighbor_elements.count('O')
            s_on_carbon = carbon_neighbor_elements.count('S')
            
            # Create detailed environment key
            if is_aromatic_carbon:
                env_key = f"aromatic_H_{carbon_hybridization}_{tuple(sorted(carbon_neighbor_elements))}"
            else:
                # Aliphatic carbon - be more specific about environment
                if h_on_carbon == 2 and c_on_carbon == 1:
                    # Check what the carbon neighbor is connected to
                    carbon_neighbor_idx = [n for n in carbon_neighbors if periodictable.elements[atomic_numbers[n]].symbol.upper() == 'C'][0]
                    carbon_neighbor_neighbors = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in neighbors[carbon_neighbor_idx]]
                    env_key = f"methyl_H_connected_to_C_with_{tuple(sorted(carbon_neighbor_neighbors))}"
                elif h_on_carbon == 1 and c_on_carbon == 2:
                    # Methylene - check what carbons are connected to
                    carbon_neighbor_envs = []
                    for cn_idx in [n for n in carbon_neighbors if periodictable.elements[atomic_numbers[n]].symbol.upper() == 'C']:
                        cn_neighbors = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in neighbors[cn_idx]]
                        carbon_neighbor_envs.append(tuple(sorted(cn_neighbors)))
                    env_key = f"methylene_H_{tuple(sorted(carbon_neighbor_envs))}"
                elif h_on_carbon == 0:
                    # Tertiary or quaternary carbon
                    env_key = f"tertiary_H_{tuple(sorted(carbon_neighbor_elements))}"
                else:
                    # General case
                    env_key = f"aliphatic_H_{carbon_hybridization}_{tuple(sorted(carbon_neighbor_elements))}"
                    
        elif heavy_element == 'N':
            # Nitrogen-bonded hydrogens
            nitrogen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            nitrogen_neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in nitrogen_neighbors]
            is_aromatic_nitrogen = heavy_atom_idx in aromatic_atoms
            
            if is_aromatic_nitrogen:
                env_key = f"aromatic_NH_{tuple(sorted(nitrogen_neighbor_elements))}"
            else:
                env_key = f"amine_H_{tuple(sorted(nitrogen_neighbor_elements))}"
                
        elif heavy_element == 'O':
            # Oxygen-bonded hydrogens
            oxygen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            oxygen_neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in oxygen_neighbors]
            env_key = f"hydroxyl_H_{tuple(sorted(oxygen_neighbor_elements))}"
            
        elif heavy_element == 'S':
            # Sulfur-bonded hydrogens
            sulfur_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            sulfur_neighbor_elements = [periodictable.elements[atomic_numbers[n]].symbol.upper() for n in sulfur_neighbors]
            env_key = f"thiol_H_{tuple(sorted(sulfur_neighbor_elements))}"
            
        else:
            # Other heavy atoms
            env_key = f"{heavy_element}_H"
        
        # Group hydrogens with identical environments
        if env_key not in hydrogen_environments:
            hydrogen_environments[env_key] = []
        hydrogen_environments[env_key].append(i)
    
    # Create constraint groups for hydrogens with identical environments
    hydrogen_constraint_groups = []
    for env_key, h_indices in hydrogen_environments.items():
        if len(h_indices) > 1:
            hydrogen_constraint_groups.append({
                'type': 'equivalent_hydrogens',
                'indices': h_indices,
                'description': f"Equivalent hydrogens: {env_key}",
                'environment': env_key
            })
    
    return hydrogen_constraint_groups

#Enhanced equivalent atom identification using chemical environment fingerprints
def identify_equivalent_atoms_enhanced(coordinates, atomic_numbers, atom_names):
    
    connectivity, neighbors, distance_matrix = analyze_molecular_connectivity(coordinates, atomic_numbers, atom_names)
    
    natom = len(atomic_numbers)
    equivalent_groups = []
    processed = set()
    
    # Group atoms by element first
    elements_dict = {}
    for i in range(natom):
        element = periodictable.elements[atomic_numbers[i]].symbol.upper()
        if element not in elements_dict:
            elements_dict[element] = []
        elements_dict[element].append(i)
    
    # For each element, find truly equivalent atoms
    for element, atom_indices in elements_dict.items():
        if len(atom_indices) < 2:
            continue
            
        # Create environment fingerprints for all atoms of this element
        fingerprints = {}
        for atom_idx in atom_indices:
            if atom_idx in processed:
                continue
                
            fingerprint = get_chemical_environment_fingerprint(atomic_numbers, atom_idx, neighbors, max_depth=3)
            
            if fingerprint not in fingerprints:
                fingerprints[fingerprint] = []
            fingerprints[fingerprint].append(atom_idx)
        
        # Create equivalent groups for atoms with identical fingerprints
        for fingerprint, equivalent_atoms in fingerprints.items():
            if len(equivalent_atoms) > 1:
                equivalent_groups.append(equivalent_atoms)
                processed.update(equivalent_atoms)
                
    # Also analyze hydrogen environments specifically
    hydrogen_groups = analyze_enhanced_hydrogen_environments(atomic_numbers, neighbors, connectivity)
    for h_group in hydrogen_groups:
        equivalent_groups.append(h_group['indices'])
    
    return equivalent_groups

#Write equivalent atoms to file in the required format
def write_equivalent_atoms(equivalent_groups, filename="similar.dat"):
    
    with open(filename, "w") as similar_file:
        for group in equivalent_groups:
            if len(group) > 1:
                # Use the first atom as reference
                reference = group[0]
                for atom in group[1:]:
                    similar_file.write(f"{atom+1:5d} {reference+1:5d}\n")  # Convert to 1-based indexing

def main(hessian_file, log_file):
    # File paths for storing computed distance, angle, and dihedral data
    distance_file = "distance.dat"
    angle_file = "angle.dat"
    dihedral_file = "dihedral.dat"
    
    # Read data from the formatted checkpoint file (fchk), extracting coordinates, Hessian, atomic numbers, and charges
    coordinates, atomic_numbers, atom_names = extract_coordinates_from_file(hessian_file)
    with open("atomic_number.dat", "w") as f:
        for number in atomic_numbers:
            f.write(f"{number}\n")

    hessian = extract_hessian_from_file(hessian_file)
    mulliken_charges = extract_mulliken_charges(log_file)
    atom_pairs = read_distances(distance_file)
    angle_definitions = read_angles(angle_file)
    dihedral_definitions = read_dihedrals(dihedral_file)
    
    # Apply the Seminario method to calculate bond force constants using Hessian, coordinates, and atomic pairs
    bonds = seminario_method(hessian, coordinates, atom_pairs, atomic_numbers)
    
    # Detect atoms with similar properties based on their coordinates, atomic numbers, and charges
    equivalent_groups = identify_equivalent_atoms_enhanced(coordinates, atomic_numbers, atom_names)
    write_equivalent_atoms(equivalent_groups)
    
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
            atom1, atom2, atom3 = angle_def[:3]
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
    if len(sys.argv) != 3:
        print("Usage: python script_name.py path/to/hessian_file path/to/log_file")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])# Run all the function
