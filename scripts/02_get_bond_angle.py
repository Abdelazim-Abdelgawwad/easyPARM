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


import numpy as np
from scipy.spatial.distance import pdist, squareform
import argparse

# Covalent radii dictionary (values in Angstroms)
covalent_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06,
    'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50,
    'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16,
    'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42,
    'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40,
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98,
    'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
    'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45,
    'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
    'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69
}

# Comman bond length to help in the detection of single, double and triple bond
typical_bond_lengths = {
    ('C', 'C'): (1.20, 1.33, 1.54),  # triple, double, single
    ('C', 'N'): (1.16, 1.27, 1.47),
    ('C', 'O'): (1.13, 1.22, 1.43),
    ('N', 'N'): (1.10, 1.25, 1.45),
    ('N', 'O'): (1.06, 1.21, 1.36),
    ('O', 'O'): (1.21, 1.48),  # double, single (no common triple bond)
    ('H', 'H'): (0.74,),  # single
    ('H', 'C'): (1.09,),  # single
    ('H', 'N'): (1.01,),  # single
    ('H', 'O'): (0.96,),  # single
    ('H', 'F'): (0.92,),  # single
    ('H', 'Cl'): (1.27,),  # single
    ('H', 'Br'): (1.41,),  # single
    ('H', 'I'): (1.61,),  # single
    ('C', 'F'): (1.35,),  # single
    ('C', 'Cl'): (1.77,),  # single
    ('C', 'Br'): (1.94,),  # single
    ('C', 'I'): (2.14,),  # single
    ('C', 'S'): (1.82,),  # single
    ('C', 'P'): (1.84,),  # single
    ('N', 'F'): (1.36,),  # single
    ('N', 'Cl'): (1.75,),  # single
    ('N', 'Br'): (1.93,),  # single
    ('N', 'I'): (2.11,),  # single
    ('N', 'S'): (1.68,),  # single
    ('N', 'P'): (1.70,),  # single
    ('O', 'F'): (1.42,),  # single
    ('O', 'Cl'): (1.70,),  # single
    ('O', 'Br'): (1.89,),  # single
    ('O', 'I'): (2.08,),  # single
    ('O', 'S'): (1.51,),  # single
    ('O', 'P'): (1.63,),  # single
    ('F', 'F'): (1.42,),  # single
    ('F', 'Cl'): (1.63,),  # single
    ('F', 'Br'): (1.76,),  # single
    ('F', 'I'): (1.91,),  # single
    ('F', 'S'): (1.56,),  # single
    ('F', 'P'): (1.54,),  # single
    ('Cl', 'Cl'): (1.99,),  # single
    ('Cl', 'Br'): (2.14,),  # single
    ('Cl', 'I'): (2.32,),  # single
    ('Cl', 'S'): (2.03,),  # single
    ('Cl', 'P'): (2.03,),  # single
    ('Br', 'Br'): (2.28,),  # single
    ('Br', 'I'): (2.47,),  # single
    ('Br', 'S'): (2.21,),  # single
    ('Br', 'P'): (2.22,),  # single
    ('I', 'I'): (2.67,),  # single
    ('I', 'S'): (2.41,),  # single
    ('I', 'P'): (2.44,),  # single
    ('S', 'S'): (2.05,),  # single
    ('S', 'P'): (2.12,),  # single
    ('P', 'P'): (2.21,),  # single
    ('O', 'P'): (1.45, 1.63),  # double, single
    ('S', 'P'): (1.95, 2.12),  # double, single
    ('N', 'P'): (1.57, 1.70),  # double, single
    ('C', 'P'): (1.69, 1.84),  # double, single
}

def is_metal(element):
    """Check if the element is a metal."""
    metals = {
        'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra',
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
        'Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi', 'Po', 'Nh', 'Fl', 'Mc', 'Lv',
        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
    }
    return element in metals

# Read the XYZ file to extract atomic coordinates and atom types.
def parse_xyz(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        atom_types = []
        coords = []
        for line in lines[2:]:  # Skip the first two lines of the XYZ file (header)
            parts = line.split()
            atom_types.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(coords), atom_types

# Detect bonded atom pairs based on covalent radii.
def detect_bonds(coords, atom_types, tolerance=0.2):
    distances = squareform(pdist(coords))
    bonds = []
    num_atoms = len(coords)
    
    def get_max_bond_length(atom_i, atom_j, base_length):
        if is_metal(atom_i) and is_metal(atom_j):
            return base_length * 1.1 
        elif is_metal(atom_i) or is_metal(atom_j):
            return base_length * (1 + tolerance )  # Increased for metal-nonmetal
        else:
            return base_length * (1 + tolerance )  # Increased for nonmetal-nonmetal

    # Calculate coordination numbers 
    def get_coordination_number(i, conservative_distances):
        return sum(1 for k in range(num_atoms) 
                  if k != i and conservative_distances[i][k])

    conservative_bonds = [[False] * num_atoms for _ in range(num_atoms)]
    for i in range(num_atoms):
        for j in range(num_atoms):
            if i != j:
                r1 = covalent_radii.get(atom_types[i], 0)
                r2 = covalent_radii.get(atom_types[j], 0)
                base_length = r1 + r2
                max_length = get_max_bond_length(atom_types[i], atom_types[j], base_length)
                conservative_bonds[i][j] = distances[i, j] <= max_length

    # Main bond detection loop
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            r1 = covalent_radii.get(atom_types[i], 0)
            r2 = covalent_radii.get(atom_types[j], 0)
            
            metal_i = is_metal(atom_types[i])
            metal_j = is_metal(atom_types[j])
            
            base_bond_length = r1 + r2
            current_distance = distances[i, j]
            
            # Determine maximum bond length
            if metal_i and metal_j:  # Metal-metal bonds
                max_bond_length = base_bond_length * 1.03
                coord_limit = 12
            elif metal_i or metal_j:  # Metal-nonmetal bonds
                max_bond_length = base_bond_length * (1 + tolerance * 1.0)  
                coord_limit = 10  # Increased limit
            else:  # Nonmetal-nonmetal bonds
                max_bond_length = base_bond_length * (1 + tolerance * 1.0)
                coord_limit = 8  # Increased limit
            
            if current_distance <= max_bond_length:
                # Get coordination numbers
                coord_num_i = get_coordination_number(i, conservative_bonds)
                coord_num_j = get_coordination_number(j, conservative_bonds)
                
                # Set coordination limits
                limit_i = coord_limit if metal_i else 8
                limit_j = coord_limit if metal_j else 8
                
                if coord_num_i <= limit_i and coord_num_j <= limit_j:
                    bonds.append((i, j))

    
    return bonds

# Detect angle triplets based on detected bonds.
def detect_angles(bonds):
    angles = []
    bond_dict = {}

    for i, j in bonds:
        if i not in bond_dict:
            bond_dict[i] = []
        if j not in bond_dict:
            bond_dict[j] = []
        bond_dict[i].append(j)
        bond_dict[j].append(i)

    for atom, neighbors in bond_dict.items():
        if len(neighbors) > 1:
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    angles.append((neighbors[i], atom, neighbors[j]))

    return angles

# Detects important dihedral based on detected bonds,
# prioritizing metal-containing and other critical dihedrals.
def detect_dihedrals(bonds, atom_types, max_dihedrals=4000):
    bond_dict = {}
    for i, j in bonds:
        if i not in bond_dict:
            bond_dict[i] = set()
        if j not in bond_dict:
            bond_dict[j] = set()
        bond_dict[i].add(j)
        bond_dict[j].add(i)

    def atom_priority(atom_idx):
        """Assign priority to atoms based on type and connectivity."""
        atom_type = atom_types[atom_idx]
        connectivity = len(bond_dict[atom_idx])
        if is_metal(atom_type):
            return 5  # Highest priority for metals
        elif atom_type in ['C', 'N', 'O', 'S', 'P']:
            return 4 if connectivity > 2 else 3
        elif atom_type != 'H':
            return 2
        return 1 if connectivity > 1 else 0  # Lower priority for terminal hydrogens

    potential_dihedrals = []
    for b in bond_dict:
        for a in bond_dict[b]:
            for c in bond_dict[b]:
                if a != c:
                    for d in bond_dict[c]:
                        if d != b:
                            priority = sum(atom_priority(atom) for atom in (a, b, c, d))
                            if any(is_metal(atom_types[atom]) for atom in (a, b, c, d)):
                                priority += 10
                            potential_dihedrals.append(((a, b, c, d), priority))

    sorted_dihedrals = sorted(set(potential_dihedrals), key=lambda x: (-x[1], x[0]))
    selected_dihedrals = [dihedral for dihedral, _ in sorted_dihedrals[:max_dihedrals]]
    
    return selected_dihedrals
    
#Calculates the distance between two atoms given their coordinates.
def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

#Calculates the angle between three atoms given their coordinates.
def calculate_angle(coord1, coord2, coord3):
    vec1 = coord1 - coord2
    vec2 = coord3 - coord2
    cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    angle = np.degrees(np.arccos(cosine_angle))
    return angle

#Calculates the dihedral angle between four atoms given their coordinates.
def calculate_dihedral(coord1, coord2, coord3, coord4):
    b1 = coord2 - coord1
    b2 = coord3 - coord2
    b3 = coord4 - coord3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    m1 = np.cross(n1, b2)

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    dihedral = np.degrees(np.arctan2(y, x))
    return dihedral

#Saves the bond distances to a file.
def save_distances(bonds, coords, filename='distance.dat'):
    with open(filename, 'w') as f:
        for bond in bonds:
            atom1, atom2 = bond
            distance = calculate_distance(coords[atom1], coords[atom2])
            f.write(f"{atom1+1}  {atom2+1}   {distance:.2f} \n")

# Saves the bond angles to a file.
def save_angles(angles, coords, filename='angle.dat'):
    with open(filename, 'w') as f:
        for angle_triplet in angles:
            atom1, atom2, atom3 = angle_triplet
            angle = calculate_angle(coords[atom1], coords[atom2], coords[atom3])
            f.write(f"{atom1+1}  {atom2+1}  {atom3+1}  {angle:.2f}\n")

#Saves the dihedral angles to a file.
def save_dihedrals(dihedrals, coords, filename='dihedral.dat'):
    with open(filename, 'w') as f:
        for dihedral_quadruplet in dihedrals:
            atom1, atom2, atom3, atom4 = dihedral_quadruplet
            angle = calculate_dihedral(coords[atom1], coords[atom2], coords[atom3], coords[atom4])
            f.write(f"{atom1+1}  {atom2+1}  {atom3+1}  {atom4+1}  {angle:.2f}\n")
    
#Saves the bond distances and types of bond to a file.
def save_distances_with_type(bonds, coords, atom_types, filename='distance_type.dat'):
    with open(filename, 'w') as f:
        for bond in bonds:
            atom1, atom2 = bond
            distance = calculate_distance(coords[atom1], coords[atom2])
            bond_type = detect_bond_type(atom1, atom2, distance, atom_types, coords, bonds)
            f.write(f"{atom1+1}  {atom2+1}  {bond_type:.1f}\n")
    
#Determine the bond type based on the elements, their distance, and molecular context.
def get_bond_type(element1, element2, distance, coords, atom_types, bonds):
    if is_metal(element1) or is_metal(element2):
        return 1

    key = tuple(sorted((element1, element2)))
    
    if key not in typical_bond_lengths:
        return 1

    lengths = typical_bond_lengths[key]
    
    tolerance_single = 1.15
    tolerance_double = 1.08
    tolerance_triple = 1.05

    if len(lengths) == 1:
        return 1
    elif len(lengths) == 2:
        if distance <= lengths[0] * tolerance_double:
            return 2
        else:
            return 1
    elif len(lengths) == 3:
        if distance <= lengths[0] * tolerance_triple:
            return 3
        elif distance <= lengths[1] * tolerance_double:
            return 2
        else:
            return 1

    return 1

#Consider the molecular context for Carbon-Carbon bonds.
def consider_carbon_context(element1, element2, distance, coords, atom_types, bonds):
    if is_part_of_aromatic_ring(element1, element2, coords, atom_types, bonds):
        return 1.5  # Represent aromatic bonds as 1.5
    if is_part_of_conjugated_system(element1, element2, coords, atom_types, bonds):
        return 1.5  # Represent conjugated bonds as 1.5 as well
    return 1

#Check if the bond is part of an aromatic ring.
def is_part_of_aromatic_ring(atom1, atom2, coords, atom_types, bonds):
    def find_ring(start, current, path, visited):
        if len(path) > 6:
            return None
        if start == current and len(path) == 6:
            return path
        for neighbor in bond_dict[current]:
            if neighbor not in visited:
                result = find_ring(start, neighbor, path + [neighbor], visited | {neighbor})
                if result:
                    return result
        return None

    bond_dict = {i: set() for i in range(len(coords))}
    for a, b in bonds:
        bond_dict[a].add(b)
        bond_dict[b].add(a)

    ring = find_ring(atom1, atom2, [atom1, atom2], {atom1, atom2})
    
    if ring:
        # Check if all atoms in the ring are carbon or nitrogen
        if all(atom_types[i] in ['C', 'N'] for i in ring):
            # Check planarity
            normal = np.cross(coords[ring[1]] - coords[ring[0]], coords[ring[2]] - coords[ring[0]])
            normal /= np.linalg.norm(normal)
            if all(abs(np.dot(normal, coords[i] - coords[ring[0]])) < 0.1 for i in ring):
                return True
    return False

#Check if the bond is part of a conjugated system.
def is_part_of_conjugated_system(atom1, atom2, coords, atom_types, bonds):
    def find_conjugated_path(start, current, path, visited):
        if len(path) > 8:  # Limit the search to avoid excessive computation
            return None
        for neighbor in bond_dict[current]:
            if neighbor not in visited:
                bond_type = get_bond_type(atom_types[current], atom_types[neighbor], 
                                          calculate_distance(coords[current], coords[neighbor]),
                                          coords, atom_types, bonds)
                if bond_type in [1.5, 2]:  # Consider aromatic (1.5) and double (2) bonds
                    result = find_conjugated_path(start, neighbor, path + [neighbor], visited | {neighbor})
                    if result:
                        return result
        return path if len(path) > 2 else None

    bond_dict = {i: set() for i in range(len(coords))}
    for a, b in bonds:
        bond_dict[a].add(b)
        bond_dict[b].add(a)

    conjugated_path = find_conjugated_path(atom1, atom2, [atom1, atom2], {atom1, atom2})
    return bool(conjugated_path)

#Get an initial guess for bond type based on element types and distance.
def get_initial_bond_type(element1, element2, distance):
    key = tuple(sorted((element1, element2)))

    if key not in typical_bond_lengths:
        return 1

    lengths = typical_bond_lengths[key]

    tolerance_single = 1.15
    tolerance_double = 1.08
    tolerance_triple = 1.05

    if len(lengths) == 1:
        return 1
    elif len(lengths) == 2:
        if distance <= lengths[0] * tolerance_double:
            return 2
        else:
            return 1
    elif len(lengths) == 3:
        if distance <= lengths[0] * tolerance_triple:
            return 3
        elif distance <= lengths[1] * tolerance_double:
            return 2
        else:
            return 1

    return 1

#Detects the bond type (single, double, or triple) based on atom types, distance, and molecular context.
def detect_bond_type(atom1, atom2, distance, atom_types, coords, bonds):
    type1, type2 = atom_types[atom1], atom_types[atom2]
    return get_bond_type(type1, type2, distance, coords, atom_types, bonds)

#Adjust bond types to ensure consistency in aromatic and conjugated systems.
def adjust_bond_types(bond_types, bonds, coords, atom_types):
    adjusted_bond_types = bond_types.copy()

    # Create a bond graph
    bond_graph = {i: set() for i in range(len(coords))}
    for bond in bonds:
        atom1, atom2 = bond
        bond_graph[atom1].add(atom2)
        bond_graph[atom2].add(atom1)

    # Detect and adjust aromatic rings
    aromatic_rings = detect_aromatic_rings(bond_graph, coords, atom_types)
    for ring in aromatic_rings:
        for i in range(len(ring)):
            atom1, atom2 = ring[i], ring[(i+1) % len(ring)]
            adjusted_bond_types[tuple(sorted((atom1, atom2)))] = 1.5

    # Adjust conjugated systems
    conjugated_systems = detect_conjugated_systems(bond_graph, bond_types, coords, atom_types)
    for system in conjugated_systems:
        adjust_conjugated_system(system, adjusted_bond_types)

    return adjusted_bond_types

#Detect aromatic rings in the molecule.
def detect_aromatic_rings(bond_graph, coords, atom_types):
    aromatic_rings = []
    visited = set()

    def dfs(atom, path):
        if len(path) > 6:
            return
        if len(path) == 6 and atom == path[0]:
            if is_aromatic(path, coords, atom_types):
                aromatic_rings.append(path)
            return
        for neighbor in bond_graph[atom]:
            if neighbor not in path[1:]:
                dfs(neighbor, path + [neighbor])

    for atom in bond_graph:
        if atom not in visited:
            dfs(atom, [atom])
            visited.add(atom)

    return aromatic_rings

#Check if a ring is aromatic.
def is_aromatic(ring, coords, atom_types):
    # Check if all atoms are C or N
    if not all(atom_types[atom] in ['C', 'N'] for atom in ring):
        return False

    # Check planarity
    normal = np.cross(coords[ring[1]] - coords[ring[0]], coords[ring[2]] - coords[ring[0]])
    normal /= np.linalg.norm(normal)
    if not all(abs(np.dot(normal, coords[atom] - coords[ring[0]])) < 0.1 for atom in ring):
        return False

    return True

#Detect conjugated systems in the molecule.
def detect_conjugated_systems(bond_graph, bond_types, coords, atom_types):
    conjugated_systems = []
    visited = set()

    def dfs(atom, system):
        visited.add(atom)
        system.append(atom)
        for neighbor in bond_graph[atom]:
            if neighbor not in visited:
                bond = tuple(sorted((atom, neighbor)))
                if bond_types[bond] in [1.5, 2]:
                    dfs(neighbor, system)

    for atom in bond_graph:
        if atom not in visited:
            system = []
            dfs(atom, system)
            if len(system) > 2:
                conjugated_systems.append(system)

    return conjugated_systems

#Adjust bond types in a conjugated system to alternate single and double bonds.
def adjust_conjugated_system(system, bond_types):
    for i in range(len(system) - 1):
        bond = tuple(sorted((system[i], system[i+1])))
        bond_types[bond] = 2 if i % 2 == 0 else 1

#Assign bond orders to all bonds in the molecule, considering aromatic systems.
def assign_bond_orders(bonds, coords, atom_types):
    bond_orders = {}
    aromatic_bonds = set()

    # First pass: Identify aromatic bonds and assign initial bond orders
    for bond in bonds:
        atom1, atom2 = bond
        distance = calculate_distance(coords[atom1], coords[atom2])
        initial_type = get_bond_type(atom_types[atom1], atom_types[atom2], distance, coords, atom_types, bonds)

        if is_part_of_aromatic_ring(atom1, atom2, coords, atom_types, bonds):
            aromatic_bonds.add(bond)
            bond_orders[bond] = 1.5
        else:
            bond_orders[bond] = initial_type

    #Second pass: Adjust bond orders to ensure no atom has more than one double bond
    for atom in range(len(coords)):
        double_bonds = [b for b in bonds if atom in b and bond_orders[b] == 2]
        if len(double_bonds) > 1:
            for bond in double_bonds[1:]:
                bond_orders[bond] = 1

    #Third pass: Assign bond orders to aromatic systems
    for bond in aromatic_bonds:
        atom1, atom2 = bond
        if all(bond_orders.get((atom1, n), 1) < 2 and bond_orders.get((atom2, n), 1) < 2
               for n in range(len(coords)) if n != atom1 and n != atom2):
            bond_orders[bond] = 1.5
        else:
            bond_orders[bond] = 1

    return bond_orders

#Saves the bond distances and types to a file.
def save_distances_with_type(bonds, coords, atom_types, filename='distance_type.dat'):
    bond_orders = assign_bond_orders(bonds, coords, atom_types)
    
    with open(filename, 'w') as f:
        for bond in bonds:
            atom1, atom2 = bond
            distance = calculate_distance(coords[atom1], coords[atom2])
            bond_type = bond_orders[bond]
            f.write(f"{atom1+1}  {atom2+1}  {bond_type}\n")

def main():
    parser = argparse.ArgumentParser(description='Process an XYZ file to calculate distances, angles, and dihedrals.')
    parser.add_argument('input_file', type=str, help='Path to the input XYZ file')
    args = parser.parse_args()

    coords, atom_types = parse_xyz(args.input_file)
    bonds = detect_bonds(coords, atom_types)
    angles = detect_angles(bonds)
    dihedrals = detect_dihedrals(bonds, atom_types)

    save_distances(bonds, coords)
    save_angles(angles, coords)
    save_dihedrals(dihedrals, coords)
    save_distances_with_type(bonds, coords, atom_types)

if __name__ == "__main__":
    main()
