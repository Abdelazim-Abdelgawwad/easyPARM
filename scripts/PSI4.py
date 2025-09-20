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

import psi4
import resp
import numpy as np
import os
import sys

bohr_to_angstrom = 0.529177249

#Read configuration parameters from psi4.config file
def read_config_file(config_file='psi4.config'):
    config = {}
    default_config = {
        'XYZ': None,
        'CHARGE': 0,
        'MULTIPLICITY': 1,
        'CPU': 4,
        'MEMORY': '4GB',
        'METHOD': 'B3LYP',
        'BASIS_SET': None  # Will use defaults if not specified
    }
    
    if not os.path.exists(config_file):
        print(f"[ERROR] Configuration file '{config_file}' not found!")
        return None
    
    try:
        with open(config_file, 'r') as f:
            lines = f.readlines()
        
        
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse key = value pairs
            if '=' in line:
                key, value = line.split('=', 1)
                key = key.strip().upper()
                value = value.strip()
                
                # Remove quotes if present
                if value.startswith('"') and value.endswith('"'):
                    value = value[1:-1]
                elif value.startswith("'") and value.endswith("'"):
                    value = value[1:-1]
                
                # Remove inline comments
                if '#' in value:
                    value = value.split('#')[0].strip()
                
                # Process different parameter types
                if key == 'XYZ':
                    config['XYZ'] = value
                elif key == 'CHARGE':
                    try:
                        config['CHARGE'] = int(value)
                    except ValueError:
                        print(f"[WARNING] Invalid CHARGE value on line {line_num}: '{value}', using default 0")
                        config['CHARGE'] = 0
                elif key == 'MULTIPLICITY':
                    try:
                        config['MULTIPLICITY'] = int(value)
                    except ValueError:
                        print(f"[WARNING] Invalid MULTIPLICITY value on line {line_num}: '{value}', using default 1")
                        config['MULTIPLICITY'] = 1
                elif key == 'CPU' or key == 'NCPUS':
                    try:
                        config['CPU'] = int(value)
                    except ValueError:
                        print(f"[WARNING] Invalid CPU value on line {line_num}: '{value}', using default 4")
                        config['CPU'] = 4
                elif key == 'MEMORY':
                    # Handle different memory formats (5GB, 500MB, etc.)
                    value_upper = value.upper()
                    if value_upper.endswith('GB') or value_upper.endswith('MB') or value_upper.endswith('KB'):
                        config['MEMORY'] = value_upper
                    else:
                        # If no unit specified, assume GB
                        try:
                            mem_val = float(value)
                            config['MEMORY'] = f"{mem_val}GB"
                        except ValueError:
                            print(f"[WARNING] Invalid MEMORY value on line {line_num}: '{value}', using default 8GB")
                            config['MEMORY'] = '8GB'
                elif key == 'METHOD':
                    config['METHOD'] = value.upper()
                elif key == 'BASIS_SET' or key == 'BASIS':
                    # Parse basis set specification
                    config['BASIS_SET'] = value
                else:
                    print(f"[WARNING] Unknown parameter on line {line_num}: '{key}'")
            else:
                print(f"[WARNING] Invalid format on line {line_num}: '{line}'")
        
        # Apply defaults for missing parameters
        for key, default_value in default_config.items():
            if key not in config:
                config[key] = default_value
                print(f"Using default {key}: {default_value}")
        
        # Validate required parameters
        if config['XYZ'] is None:
            print("[ERROR] XYZ file not specified in configuration!")
            return None
        
        # Check if XYZ file exists
        if not os.path.exists(config['XYZ']):
            print(f"[ERROR] XYZ file '{config['XYZ']}' not found!")
            return None
        
        if config['BASIS_SET']:
            print(f"  Basis set: {config['BASIS_SET']}")
        else:
            print(f"  Basis set: Default (6-31G* for non-metals, LANL2DZ for metals)")
        
        return config
        
    except Exception as e:
        print(f"[ERROR] Failed to read configuration file: {e}")
        return None

#Read XYZ file and return coordinates as string
def read_xyz_file(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"XYZ file '{filename}' not found!")
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Skip the first two lines (number of atoms and comment)
        if len(lines) < 3:
            raise ValueError(f"Invalid XYZ file format: {filename}")
        
        try:
            natoms = int(lines[0].strip())
        except ValueError:
            raise ValueError(f"First line of XYZ file must contain number of atoms")
        
        if len(lines) < natoms + 2:
            raise ValueError(f"XYZ file has insufficient lines for {natoms} atoms")
        
        # Extract coordinate lines
        coord_lines = []
        for i in range(2, min(len(lines), natoms + 2)):
            line = lines[i].strip()
            if line:  # Skip empty lines
                coord_lines.append(line)
        
        xyz_coords = '\n'.join(coord_lines)
        
        return xyz_coords
        
    except Exception as e:
        raise Exception(f"Failed to read XYZ file '{filename}': {e}")

#Parse XYZ coordinates and clean them up for Psi4
def parse_xyz_coordinates(xyz_coordinates):
    lines = xyz_coordinates.strip().split('\n')
    clean_lines = []
    unique_elements = set()
    element_count = {}
     
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        parts = line.split()
        if len(parts) < 4:
            print(f"Warning: Line {i+1} has insufficient data: '{line}'")
            continue
            
        try:
            element = parts[0].strip()
            x = float(parts[1])
            y = float(parts[2])  
            z = float(parts[3])
            
            element = element[0].upper() + element[1:].lower() if len(element) > 1 else element.upper()
            
            unique_elements.add(element)
            element_count[element] = element_count.get(element, 0) + 1
            
            clean_lines.append(f"{element:2s} {x:15.10f} {y:15.10f} {z:15.10f}")
            
        except (ValueError, IndexError) as e:
            print(f"Error parsing line {i+1}: '{line}' - {e}")
            continue
    
    if not clean_lines:
        raise ValueError("No valid coordinates found in input")
    
    clean_coords = "\n".join(clean_lines)
    
    
    return clean_coords, unique_elements, element_count


#Save optimized geometry in proper XYZ format
def save_optimized_geometry(wfn, filename='optimized.xyz'):
    try:
        opt_mol = wfn.molecule()
        natom = opt_mol.natom()
        xyz_lines = [f"{natom}", f"Optimized geometry : energy = {wfn.energy():.10f} Hartree"]
        for i in range(natom):
            sym = opt_mol.symbol(i).capitalize()
            x, y, z = opt_mol.x(i)*bohr_to_angstrom, opt_mol.y(i)*bohr_to_angstrom, opt_mol.z(i)*bohr_to_angstrom
            xyz_lines.append(f"{sym:2} {x:15.10f} {y:15.10f} {z:15.10f}")
        content = "\n".join(xyz_lines) + "\n"
        with open(filename, 'w') as f:
            f.write(content)
        return content
    except Exception as e:
        print(f"[ERROR] Failed to save optimized geometry: {e}")
        return None

#Save Hessian matrix to text file
def save_hessian_matrix(hessian, filename='hessian.txt'):
    try:
        # Convert to numpy array if it isn't already
        if hasattr(hessian, 'to_array'):
            hess_array = hessian.to_array()
        else:
            hess_array = np.array(hessian)
        
        # Save with header information
        with open(filename, 'w') as f:
            
            # Save in a readable format with scientific notation
            for i in range(hess_array.shape[0]):
                for j in range(hess_array.shape[1]):
                    f.write(f"{hess_array[i,j]:20.12E}")
                    if j < hess_array.shape[1] - 1:
                        f.write("  ")
                f.write("\n")
        
        return hess_array
        
    except Exception as e:
        print(f"[ERROR] Failed to save Hessian matrix: {e}")
        return None

#Analyze molecular connectivity to identify chemically equivalent atoms
#Returns a dictionary mapping atom types to lists of equivalent atom indices
def analyze_molecular_connectivity(mol):
     
    natom = mol.natom()
    
    # Calculate distance matrix
    distance_matrix = np.zeros((natom, natom))
    for i in range(natom):
        for j in range(natom):
            if i != j:
                dx = mol.x(i) - mol.x(j)
                dy = mol.y(i) - mol.y(j) 
                dz = mol.z(i) - mol.z(j)
                distance_matrix[i][j] = np.sqrt(dx*dx + dy*dy + dz*dz) * bohr_to_angstrom
    
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
        
        elem1_is_metal = elem1.upper() in metal_atoms
        elem2_is_metal = elem2.upper() in metal_atoms
        
        if elem1_is_metal or elem2_is_metal:
            if elem1_is_metal and elem2_is_metal:
                return bond_cutoffs.get(('Metal', 'Metal'), 3.0)
            else:
                non_metal = elem2 if elem1_is_metal else elem1
                return bond_cutoffs.get(('Metal', non_metal.upper()), 
                                      bond_cutoffs.get((non_metal.upper(), 'Metal'), 2.5))
        else:
            pair1 = (elem1.upper(), elem2.upper())
            pair2 = (elem2.upper(), elem1.upper())
            return bond_cutoffs.get(pair1, bond_cutoffs.get(pair2, 1.8))
    
    # Build connectivity matrix
    connectivity = np.zeros((natom, natom), dtype=bool)
    for i in range(natom):
        for j in range(i+1, natom):
            elem_i = mol.symbol(i)
            elem_j = mol.symbol(j)
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

#Create a detailed chemical fingerprint for an atom based on its environment like which atom and neighbour are linked
def get_chemical_environment_fingerprint(mol, atom_idx, neighbors, max_depth=3):
    
    def explore_environment(current_atom, depth, visited, path):
        if depth > max_depth:
            return []
        
        visited.add(current_atom)
        environment_info = []
        
        # Get immediate neighbors
        current_neighbors = []
        for neighbor_idx in neighbors[current_atom]:
            if neighbor_idx not in visited:
                neighbor_element = mol.symbol(neighbor_idx).upper()
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
    fingerprint = explore_environment(atom_idx, 0, set(), [mol.symbol(atom_idx).upper()])
    return tuple(sorted(fingerprint))

#Identify aromatic rings in the molecule using a simple ring detection and planarity/conjugation heuristics.
def identify_aromatic_rings(mol, connectivity, neighbors):
    
    natom = mol.natom()
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
            carbon_count = sum(1 for atom in ring if mol.symbol(atom).upper() == 'C')
            nitrogen_count = sum(1 for atom in ring if mol.symbol(atom).upper() == 'N')
            
            # Heuristic: if mostly C and N, likely aromatic
            if (carbon_count + nitrogen_count) >= len(ring) * 0.8:
                aromatic_atoms.update(ring)
    
    return aromatic_atoms

#Estimate carbon hybridization based on bonding pattern.
def analyze_carbon_hybridization(mol, carbon_idx, neighbors):
    
    carbon_neighbors = neighbors[carbon_idx]
    num_neighbors = len(carbon_neighbors)
    
    # Count different atom types
    neighbor_elements = [mol.symbol(n).upper() for n in carbon_neighbors]
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

#Hydrogen environment analysis that considers all the envi
def analyze_enhanced_hydrogen_environments(mol, neighbors, connectivity):
    
    natom = mol.natom()
    
    # First identify aromatic rings
    aromatic_atoms = identify_aromatic_rings(mol, connectivity, neighbors)
    
    # Categorize hydrogens with detailed environment analysis
    hydrogen_environments = {}
    
    for i in range(natom):
        if mol.symbol(i).upper() != 'H':
            continue
        
        # Get the heavy atom this H is bonded to
        heavy_neighbors = []
        for neighbor_idx in neighbors[i]:
            neighbor_element = mol.symbol(neighbor_idx).upper()
            if neighbor_element != 'H':
                heavy_neighbors.append((neighbor_idx, neighbor_element))
        
        if len(heavy_neighbors) != 1:
            continue  # Skip unusual cases
            
        heavy_atom_idx, heavy_element = heavy_neighbors[0]
        
        # Create detailed environment fingerprint
        env_fingerprint = get_chemical_environment_fingerprint(mol, i, neighbors, max_depth=3)
        
        # Analyze specific environments
        if heavy_element == 'C':
            carbon_hybridization = analyze_carbon_hybridization(mol, heavy_atom_idx, neighbors)
            is_aromatic_carbon = heavy_atom_idx in aromatic_atoms
            
            # Get carbon's neighbors (excluding this hydrogen)
            carbon_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            carbon_neighbor_elements = [mol.symbol(n).upper() for n in carbon_neighbors]
            
            # Count different types of neighbors
            h_on_carbon = carbon_neighbor_elements.count('H')
            c_on_carbon = carbon_neighbor_elements.count('C')
            n_on_carbon = carbon_neighbor_elements.count('N')
            o_on_carbon = carbon_neighbor_elements.count('O')
            s_on_carbon = carbon_neighbor_elements.count('S')
            
            # Create detailed environment key
            if is_aromatic_carbon:
                env_key = f"aromatic_H_{carbon_hybridization}_{tuple(sorted(carbon_neighbor_elements))}_{env_fingerprint}"
            else:
                # Aliphatic carbon - be more specific about environment
                if h_on_carbon == 2 and c_on_carbon == 1:
                    # Check what the carbon neighbor is connected to
                    carbon_neighbor_idx = [n for n in carbon_neighbors if mol.symbol(n).upper() == 'C'][0]
                    carbon_neighbor_neighbors = [mol.symbol(n).upper() for n in neighbors[carbon_neighbor_idx]]
                    env_key = f"methyl_H_connected_to_C_with_{tuple(sorted(carbon_neighbor_neighbors))}"
                elif h_on_carbon == 1 and c_on_carbon == 2:
                    # Methylene - check what carbons are connected to
                    carbon_neighbor_envs = []
                    for cn_idx in [n for n in carbon_neighbors if mol.symbol(n).upper() == 'C']:
                        cn_neighbors = [mol.symbol(n).upper() for n in neighbors[cn_idx]]
                        carbon_neighbor_envs.append(tuple(sorted(cn_neighbors)))
                    env_key = f"methylene_H_{tuple(sorted(carbon_neighbor_envs))}"
                elif h_on_carbon == 0:
                    # Tertiary or quaternary carbon
                    env_key = f"tertiary_H_{tuple(sorted(carbon_neighbor_elements))}_{env_fingerprint}"
                else:
                    # General case
                    env_key = f"aliphatic_H_{carbon_hybridization}_{tuple(sorted(carbon_neighbor_elements))}_{env_fingerprint}"
                    
        elif heavy_element == 'N':
            # Nitrogen-bonded hydrogens
            nitrogen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            nitrogen_neighbor_elements = [mol.symbol(n).upper() for n in nitrogen_neighbors]
            is_aromatic_nitrogen = heavy_atom_idx in aromatic_atoms
            
            if is_aromatic_nitrogen:
                env_key = f"aromatic_NH_{tuple(sorted(nitrogen_neighbor_elements))}_{env_fingerprint}"
            else:
                env_key = f"amine_H_{tuple(sorted(nitrogen_neighbor_elements))}_{env_fingerprint}"
                
        elif heavy_element == 'O':
            # Oxygen-bonded hydrogens
            oxygen_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            oxygen_neighbor_elements = [mol.symbol(n).upper() for n in oxygen_neighbors]
            env_key = f"hydroxyl_H_{tuple(sorted(oxygen_neighbor_elements))}_{env_fingerprint}"
            
        elif heavy_element == 'S':
            # Sulfur-bonded hydrogens
            sulfur_neighbors = [n for n in neighbors[heavy_atom_idx] if n != i]
            sulfur_neighbor_elements = [mol.symbol(n).upper() for n in sulfur_neighbors]
            env_key = f"thiol_H_{tuple(sorted(sulfur_neighbor_elements))}_{env_fingerprint}"
            
        else:
            # Other heavy atoms
            env_key = f"{heavy_element}_H_{env_fingerprint}"
        
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

#Equivalent atom identification that uses chemical environment fingerprints
def identify_equivalent_atoms_enhanced(mol, connectivity, neighbors):
    
    natom = mol.natom()
    equivalent_groups = []
    processed = set()
    
    # Group atoms by element first
    elements_dict = {}
    for i in range(natom):
        element = mol.symbol(i).upper()
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
                
            fingerprint = get_chemical_environment_fingerprint(mol, atom_idx, neighbors, max_depth=3)
            
            if fingerprint not in fingerprints:
                fingerprints[fingerprint] = []
            fingerprints[fingerprint].append(atom_idx)
        
        # Create equivalent groups for atoms with identical fingerprints
        for fingerprint, equivalent_atoms in fingerprints.items():
            if len(equivalent_atoms) > 1:
                equivalent_groups.append(equivalent_atoms)
                processed.update(equivalent_atoms)
                
                # Print detailed information about the group
                elements = [mol.symbol(idx).capitalize() for idx in equivalent_atoms]
    
    return equivalent_groups

#Identify special hydrogen environments that should be constrained together
def analyze_special_hydrogen_environments(mol, neighbors):
    natom = mol.natom()
    special_h_groups = []
    
    # Find different types of hydrogen environments
    h_environments = {
        'aromatic_h': [],      # H bonded to aromatic carbons
        'methyl_h': [],        # H in methyl groups
        'methylene_h': [],     # H in methylene groups
        'amine_h': [],         # H bonded to nitrogen
        'hydroxyl_h': [],      # H bonded to oxygen
        'thiol_h': []          # H bonded to sulfur
    }
    
    for i in range(natom):
        if mol.symbol(i).upper() != 'H':
            continue
            
        # Get the heavy atom this H is bonded to
        heavy_neighbors = []
        for neighbor_idx in neighbors[i]:
            neighbor_element = mol.symbol(neighbor_idx).upper()
            if neighbor_element != 'H':
                heavy_neighbors.append((neighbor_idx, neighbor_element))
        
        if len(heavy_neighbors) != 1:
            continue  # Skip if H has multiple heavy atom neighbors or no neighbors
            
        heavy_atom_idx, heavy_element = heavy_neighbors[0]
        
        # Classify hydrogen environment
        if heavy_element == 'C':
            # Count other neighbors of the carbon
            carbon_neighbors = []
            for c_neighbor in neighbors[heavy_atom_idx]:
                if c_neighbor != i:  # Exclude this hydrogen
                    carbon_neighbors.append(mol.symbol(c_neighbor).upper())
            
            h_count_on_carbon = carbon_neighbors.count('H')
            c_count_on_carbon = carbon_neighbors.count('C')
            
            if h_count_on_carbon == 2 and c_count_on_carbon >= 1:
                h_environments['methyl_h'].append(i)
            elif h_count_on_carbon == 1:
                h_environments['methylene_h'].append(i)
            else:
                # Could be aromatic or other - simple check
                if c_count_on_carbon >= 2:
                    h_environments['aromatic_h'].append(i)
                else:
                    h_environments['methyl_h'].append(i)
                    
        elif heavy_element == 'N':
            h_environments['amine_h'].append(i)
        elif heavy_element == 'O':
            h_environments['hydroxyl_h'].append(i)
        elif heavy_element == 'S':
            h_environments['thiol_h'].append(i)
    
    # Create groups for hydrogens in same environment
    for env_type, h_indices in h_environments.items():
        if len(h_indices) > 1:
            special_h_groups.append({
                'type': env_type,
                'indices': h_indices,
                'description': f"{env_type.replace('_', ' ').title()} hydrogens"
            })
    
    return special_h_groups

#Calculate RESP charges using the resp module with automatic constraint detection

def create_resp_constraints_enhanced(mol, equivalent_groups, special_h_groups, first_stage_charges):
    
    all_constraint_groups = []
    
    # Add general equivalent groups (non-hydrogen atoms mostly)
    for group in equivalent_groups:
        if len(group) > 1:
            # Check if group contains only non-hydrogen atoms or mixed
            elements = [mol.symbol(idx).upper() for idx in group]
            if 'H' not in elements:  # Non-hydrogen equivalent groups
                all_constraint_groups.append(group)
            else:
                # For groups containing hydrogens, we'll handle them separately
                # with the enhanced hydrogen analysis
                non_h_atoms = [idx for idx in group if mol.symbol(idx).upper() != 'H']
                if len(non_h_atoms) > 1:
                    all_constraint_groups.append(non_h_atoms)
    
    # Add hydrogen constraint groups from enhanced analysis
    already_constrained = set()
    for group in all_constraint_groups:
        already_constrained.update(group)
    
    for h_group in special_h_groups:
        h_indices = h_group['indices']
        # Only add if not already constrained
        new_indices = [idx for idx in h_indices if idx not in already_constrained]
        if len(new_indices) > 1:
            all_constraint_groups.append(new_indices)
            already_constrained.update(new_indices)
    
    # Create RESP constraint format
    constraint_groups_1based = []
    constraint_charge = []
    
    # First, collect all constrained atoms
    all_constrained_atoms = set()
    for group in all_constraint_groups:
        if len(group) >= 2:
            all_constrained_atoms.update(group)
    
    # Create constraint groups (1-based indexing for RESP)
    for group_idx, atom_indices in enumerate(all_constraint_groups):
        if len(atom_indices) < 2:
            continue
        
        # Convert to 1-based indexing for RESP
        indices_1based = [idx + 1 for idx in atom_indices]
        constraint_groups_1based.append(indices_1based)
        
        # Print group info
        elements = [mol.symbol(idx).capitalize() for idx in atom_indices]
        avg_charge = np.mean([first_stage_charges[idx] for idx in atom_indices])
    
    # Create constraint_charge list for NON-CONSTRAINED atoms only
    # These are the "free" atoms that are not part of any equivalent group
    natom = mol.natom()
    for atom_idx in range(natom):
        if atom_idx not in all_constrained_atoms:
            # This atom is not constrained, so it gets its own charge constraint
            constraint_charge.append([first_stage_charges[atom_idx], [atom_idx + 1]])  # 1-based indexing
    
    for i, (charge_val, atom_list) in enumerate(constraint_charge):
        atom_idx_0based = atom_list[0] - 1  # Convert back to 0-based for display
        element = mol.symbol(atom_idx_0based).capitalize()
    
    return constraint_groups_1based, constraint_charge

#RESP charges calculation with improved constraint detection
def calculate_resp_charges_integrated_enhanced(wfn, filename='resp_charges.dat', config_basis=None):
    
    mol = wfn.molecule()
    natom = mol.natom()
    
    metal_vdw_radii = {
        'AR': 1.88, 'K': 1.36, 'CA': 1.31, 'SC': 1.29, 'TI': 1.23,
        'V': 1.21, 'CR': 1.23, 'MN': 1.23, 'FE': 1.22, 'CO': 1.20,
        'NI': 1.20, 'CU': 1.19, 'ZN': 1.20, 'GA': 1.16, 'GE': 2.11,
        'AS': 1.85, 'SE': 1.90, 'BR': 1.85, 'KR': 2.02,
        'RB': 1.61, 'SR': 1.42, 'Y': 1.37, 'ZR': 1.26, 'NB': 1.28,
        'MO': 1.23, 'TC': 1.22, 'RU': 1.23, 'RH': 1.22, 'PD': 1.08,
        'AG': 1.26, 'CD': 1.22, 'IN': 1.22, 'SN': 1.21, 'SB': 2.06,
        'TE': 2.06, 'I': 1.98, 'XE': 2.16, 'CS': 1.74, 'BA': 1.52,
        'LA': 1.49, 'CE': 1.47, 'PR': 1.46, 'ND': 1.46, 'PM': 2.43,
        'SM': 2.42, 'EU': 2.40, 'GD': 2.38, 'TB': 2.37, 'DY': 2.35,
        'HO': 2.33, 'ER': 2.32, 'TM': 2.30, 'YB': 2.28, 'LU': 2.27,
        'HF': 1.32, 'TA': 1.27, 'W': 1.29, 'RE': 1.25, 'OS': 1.24,
        'IR': 1.21, 'PT': 1.16, 'AU': 1.16, 'HG': 1.23, 'TL': 1.24,
        'PB': 1.30, 'BI': 2.07, 'PO': 1.97, 'AT': 2.02, 'RN': 2.20,
        'FR': 3.48, 'RA': 2.83, 'AC': 2.60, 'TH': 2.37, 'PA': 2.43,
        'U': 2.40, 'NP': 2.39, 'PU': 2.43, 'AM': 2.44, 'CM': 2.45,
        'BK': 2.46, 'CF': 2.48, 'ES': 2.49, 'FM': 2.51, 'MD': 2.53,
        'NO': 2.54, 'LR': 2.56, 'RF': 2.57, 'DB': 2.59, 'SG': 2.60,
        'BH': 2.62, 'HS': 2.63, 'MT': 2.64, 'DS': 2.65, 'RG': 2.66,
        'CN': 2.68, 'NH': 2.69, 'FL': 2.70, 'MC': 2.71, 'LV': 2.72,
        'TS': 2.73, 'OG': 2.74
    }

    metal_atoms = {
        'RU', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MN', 'CR', 'V', 'TI', 'SC',
        'PD', 'PT', 'AU', 'AG', 'CD', 'HG', 'PB', 'SN', 'AL', 'GA', 'IN', 'TL',
        'MG', 'CA', 'SR', 'BA', 'LI', 'NA', 'K', 'RB', 'CS', 'Y', 'ZR', 'NB', 
        'MO', 'TC', 'RH', 'IR', 'OS', 'RE', 'W', 'TA', 'HF', 'LA', 'CE', 'PR', 'ND'
    }
    
    # Get unique elements in the molecule
    unique_elements = set()
    for i in range(mol.natom()):
        element = mol.symbol(i)
        unique_elements.add(element.upper())
    
    # Check if metals are present
    metals_found = bool(unique_elements.intersection(metal_atoms))
    
    # Determine basis set for ESP calculation (same logic as before)
    if config_basis:
        basis_spec = parse_basis_set_specification(config_basis)
        if basis_spec and metals_found:
            non_metal_basis = basis_spec['non_metal']
            metal_basis = basis_spec['metal']
            
            basis_assignments = []
            for element in unique_elements:
                if element in metal_atoms:
                    basis_assignments.append(f"assign {element} {metal_basis}")
                else:
                    basis_assignments.append(f"assign {element} {non_metal_basis}")
            
            basis_block = "\n".join(basis_assignments)
            
            try:
                psi4.basis_helper(basis_block, name="resp_user_basis")
                basis_esp = 'resp_user_basis'
            except Exception as e:
                basis_esp = '6-31G'
        elif basis_spec and not metals_found:
            basis_esp = basis_spec['non_metal']
        else:
            if metals_found:
                basis_esp = '6-31G'
            else:
                basis_esp = '6-31G*'
    else:
        # Use default logic for basis selection
        if metals_found:
            basis_assignments = []
            for element in unique_elements:
                if element in metal_atoms:
                    basis_assignments.append(f"assign {element} LANL2DZ")
                else:
                    basis_assignments.append(f"assign {element} 6-31G*")
            
            basis_block = "\n".join(basis_assignments)
            for assignment in basis_assignments:
                print(f"  {assignment}")
            try:
                psi4.basis_helper(basis_block, name="resp_mixed_basis")
                basis_esp = 'resp_mixed_basis'
            except Exception as e:
                basis_esp = '6-31G'
        else:
            basis_esp = '6-31G*'
    
    # RESP options
    options = {
        'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0],
        'VDW_POINT_DENSITY': 1.0,
        'RESP_A': 0.0005,
        'RESP_B': 0.1,
        'BASIS_ESP': basis_esp,
        'VDW_RADII': metal_vdw_radii,
    }
    
    try:
        print("Running first stage RESP fit...")
        # Call for first stage fit
        charges1 = resp.resp([mol], options)
        
        print('First stage - Electrostatic Potential Charges:')
        print(charges1[0])
        print('First stage - Restrained Electrostatic Potential Charges:')
        print(charges1[1])
        
        # Enhanced connectivity and environment analysis
        connectivity, neighbors, distance_matrix = analyze_molecular_connectivity(mol)
        
        # Use enhanced equivalent atom identification
        equivalent_groups = identify_equivalent_atoms_enhanced(mol, connectivity, neighbors)
        
        # Use enhanced hydrogen environment analysis
        special_h_groups = analyze_enhanced_hydrogen_environments(mol, neighbors, connectivity)
        
        # Enhanced constraint creation
        constraint_groups_1based, constraint_charge = create_resp_constraints_enhanced(
            mol, equivalent_groups, special_h_groups, charges1[1]
        )
        # SET UP SECOND STAGE WITH ENHANCED CONSTRAINTS
        print("\nRunning second stage RESP fit with constraints...")
        
        # Change the value of the RESP parameter A for second stage
        options['RESP_A'] = 0.001
        
        if constraint_groups_1based:
            options["constraint_charge"] = constraint_charge
            options["constraint_group"] = constraint_groups_1based
        else:
            print("No constraint groups found - running second stage without constraints")
        
        # Set up grid reuse
        options['grid'] = [f'1_{mol.name()}_grid.dat']
        options['esp'] = [f'1_{mol.name()}_grid_esp.dat']
        
        # Call for second stage fit
        charges2 = resp.resp([mol], options)
        
        print("\nSecond stage - Final RESP Charges:")
        final_charges = charges2[1]
        print(final_charges)
        
        # Save charges to file with detailed information
        with open(filename, 'w') as f:
            
            total_charge = 0.0
            for i in range(natom):
                element = mol.symbol(i).capitalize()
                charge = final_charges[i]
                total_charge += charge
                
                # Identify which constraint group this atom belongs to (if any)
                group_info = "None"
                for group_idx, group in enumerate(constraint_groups_1based):
                    if (i + 1) in group:  # Convert to 1-based for comparison
                        group_info = f"Group_{group_idx+1}"
                        break
                
                f.write(f"{i+1:5d} {element:>8} {charge:12.6f} \n")
            
            # Write constraint group details
            if constraint_groups_1based:
                for group_idx, group in enumerate(constraint_groups_1based):
                    elements_in_group = [mol.symbol(idx-1).capitalize() for idx in group]  # Convert back to 0-based
                    avg_charge = np.mean([final_charges[idx-1] for idx in group])
        
        return final_charges, "Two_stage_RESP_with_enhanced_constraints"
        
    except Exception as e:
        print(f"[ERROR] Enhanced RESP calculation failed: {e}")
        print("Falling back to Mulliken charges...")
        
        # Fallback to Mulliken charges (same as before)
        try:
            psi4.oeprop(wfn, 'MULLIKEN_CHARGES')
            mulliken_charges = wfn.atomic_point_charges()
            
            if mulliken_charges is not None:
                if hasattr(mulliken_charges, '__len__') and len(mulliken_charges) > 0:
                    charges = [float(charge) for charge in mulliken_charges]
                    
                    # Save Mulliken charges
                    with open(filename, 'w') as f:
                        f.write("# Mulliken Charges (RESP calculation failed)\n")
                        f.write("# Atom  Element   Charge\n")
                        f.write("#----  --------  ---------\n")
                        
                        for i in range(natom):
                            element = mol.symbol(i)
                            charge = charges[i]
                            f.write(f"{i+1:5d} {element:>8} {charge:15.6f}\n")
                    
                    return charges, "Mulliken_charges"
                else:
                    print("[ERROR] Mulliken charges array is empty")
                    return None, None
            else:
                print("[ERROR] Could not calculate any charges")
                return None, None
                
        except Exception as fallback_error:
            print(f"[ERROR] Mulliken charges also failed: {fallback_error}")
            return None, None

#Parse basis set specification from config file
def parse_basis_set_specification(basis_spec):
    if not basis_spec:
        return None
    
    # Split by comma and clean up
    basis_parts = [b.strip().upper() for b in basis_spec.split(',')]
    
    basis_dict = {}
    
    if len(basis_parts) == 1:
        # Single basis set for all atoms
        basis_dict['non_metal'] = basis_parts[0]
        basis_dict['metal'] = basis_parts[0]
    elif len(basis_parts) == 2:
        # Two basis sets: non-metal, metal
        basis_dict['non_metal'] = basis_parts[0]
        basis_dict['metal'] = basis_parts[1]
    elif len(basis_parts) >= 3:
        # Three or more: non-metal, metal, fallback
        basis_dict['non_metal'] = basis_parts[0]
        basis_dict['metal'] = basis_parts[1]
        basis_dict['fallback'] = basis_parts[2]
    
    return basis_dict

#Set up basis sets for calculation based on user specification or defaults
def setup_basis_sets(unique_elements, metal_atoms, config_basis=None):
    
    # Parse user-specified basis sets
    if config_basis:
        basis_spec = parse_basis_set_specification(config_basis)
        if basis_spec:
            non_metal_basis = basis_spec['non_metal']
            metal_basis = basis_spec['metal']
            fallback_basis = basis_spec.get('fallback', 'DEF2-SVP')
        else:
            # Default if parsing failed
            non_metal_basis = '6-31G*'
            metal_basis = 'LANL2DZ'
            fallback_basis = 'DEF2-SVP'
    else:
        # Use defaults
        non_metal_basis = '6-31G*'
        metal_basis = 'LANL2DZ'
        fallback_basis = 'DEF2-SVP'
    
    # Check if we need mixed basis sets
    metals_found = [elem for elem in unique_elements if elem in metal_atoms]
    non_metals_found = [elem for elem in unique_elements if elem not in metal_atoms]
    
    if metals_found and non_metals_found:
        # Mixed system - create basis assignments
 
        basis_assignments = []
        for element in unique_elements:
            if element in metal_atoms:
                basis_assignments.append(f"assign {element} {metal_basis}")
            else:
                basis_assignments.append(f"assign {element} {non_metal_basis}")
        
        basis_block = "\n".join(basis_assignments)
        basis_name = "mixed_basis"
        
        try:
            psi4.basis_helper(basis_block, name=basis_name)
            return basis_name, basis_block, basis_assignments
            
        except Exception as e:
            return fallback_basis, None, []
    
    elif metals_found and not non_metals_found:
        # Only metals
        return metal_basis, None, []
    
    elif non_metals_found and not metals_found:
        # Only non-metals
        return non_metal_basis, None, []
    
    else:
        # Fallback case
        print(f"\nUsing fallback basis set: {fallback_basis}")
        return fallback_basis, None, []
    
#Complete workflow: Optimize -> Calculate Hessian -> Calculate RESP charges
def optimize_metal_complex_with_resp(xyz_coordinates, charge=0, multiplicity=1, 
                                   memory="8GB", ncpus=6, method="B3LYP", basis_spec=None):
    
    # Set computational resources
    psi4.set_memory(memory)
    psi4.set_num_threads(ncpus)
    
    # Define metal atoms
    metal_atoms = {
        'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
        'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
    }
    
    # Parse coordinates
    try:
        clean_coords, unique_elements, element_count = parse_xyz_coordinates(xyz_coordinates)
    except Exception as e:
        print(f"Error parsing coordinates: {e}")
        return None, None, None, None
    
    metals_found = []
    organics_found = []
    
    for element in sorted(unique_elements):
        count = element_count[element]
        if element in metal_atoms:
            metals_found.append(element)
        else:
            organics_found.append(element)
    
    # Create molecule
    mol_string = f"""
{charge} {multiplicity}
{clean_coords}
units angstrom
symmetry c1
"""
    
    try:
        mol = psi4.geometry(mol_string)
        print("Molecule object created successfully")
    except Exception as e:
        print(f"Failed to create molecule object: {e}")
        return None, None, None, None
    
    # Set up basis sets using the new system
    basis_name, basis_block, basis_assignments = setup_basis_sets(
        unique_elements, metal_atoms, basis_spec
    )
    
    # Apply the basis set
    psi4.set_options({'basis': basis_name})
    
    if basis_assignments:
        for assignment in basis_assignments:
            print(f"  {assignment}")
    
    # Set calculation options
    psi4.set_options({
        'reference': 'UHF' if multiplicity > 1 else 'RHF',
        'scf_type': 'DF',
        'maxiter': 300,
        'guess': 'SAD',
        'e_convergence': 1e-8,
        'd_convergence': 1e-7,
        'g_convergence': 'GAU',
        'geom_maxiter': 300,
        'dft_nuclear_scheme': 'TREUTLER'
    })
     
    try:
        # Perform geometry optimization
        energy, wfn = psi4.optimize(method, molecule=mol, return_wfn=True)
        
        # Save optimized geometry - this is the crucial part that was failing
        optimized_geom = save_optimized_geometry(wfn, 'optimized.xyz')
        
        if optimized_geom is None:
            print("[ERROR] Failed to save optimized geometry")
            # Try alternative method
            try:
                opt_mol = wfn.molecule()
                natom = opt_mol.natom()
                with open('optimized.xyz', 'w') as f:
                    f.write(f"{natom}\n")
                    f.write(f"Optimized geometry : energy = {energy:.10f} Hartree\n")
                    for i in range(natom):
                        sym = opt_mol.symbol(i)
                        x = opt_mol.x(i) * bohr_to_angstrom
                        y = opt_mol.y(i) * bohr_to_angstrom
                        z = opt_mol.z(i) * bohr_to_angstrom
                        f.write(f"{sym:2} {x:15.10f} {y:15.10f} {z:15.10f}\n")
                optimized_geom = "Alternative_save_successful"
            except Exception as alt_error:
                print(f"[ERROR] Alternative geometry save also failed: {alt_error}")
                optimized_geom = None
        
        # Get dipole moment
        try:
            dipole = psi4.variable('CURRENT DIPOLE')
            if hasattr(dipole, '__len__') and len(dipole) > 0:
                dipole_magnitude = np.sqrt(sum(d**2 for d in dipole))
        except Exception as dipole_error:
            print(f"Could not retrieve dipole moment: {dipole_error}")
        
        hess_array = None
        try:
            hessian = psi4.hessian(method, molecule=wfn.molecule())
            
            hess_array = save_hessian_matrix(hessian, 'hessian.txt')
            
            if hess_array is not None:
                print(f"  Matrix dimensions: {hess_array.shape}")
                
        except Exception as hess_error:
            print(f"[ERROR] Hessian calculation failed: {hess_error}")
        
        # Calculate RESP charges with proper error handling
        charges = None
        method_used = None
        try:
            charges, method_used = calculate_resp_charges_integrated_enhanced(wfn, 'resp_charges.dat', basis_spec)
        except Exception as resp_error:
            print(f"[ERROR] Charge calculation failed completely: {resp_error}")
        
        if optimized_geom is not None:
            print(f"Optimized geometry saved to: optimized.xyz")
        else:
            print(f"Failed to save optimized geometry")
        
        if hess_array is not None:
            print(f"Hessian calculation: SUCCESSFUL")
            print(f"Hessian matrix saved to: hessian.txt")
        else:
            print(f"Hessian calculation: FAILED")
        
        if charges is not None and method_used is not None:
            print(f"Charges calculated using: {method_used}")
            print(f"Charges saved to: resp_charges.dat")
        else:
            print(f"Charge calculation: FAILED")
        
        return energy, optimized_geom, hess_array, charges
        
    except Exception as e:
        
        # Print more detailed error information
        import traceback
        print("Detailed traceback:")
        traceback.print_exc()
        
        return None, None, None, None

#Create a sample configuration file for reference
def create_sample_config(filename='psi4.config'):
    sample_content = """# Psi4 Configuration File
# Lines starting with # are comments
# Format: PARAMETER = VALUE

# XYZ file containing the molecular structure (required)
XYZ = structure.xyz

# Molecular charge (integer)
CHARGE = 0

# Spin multiplicity 
MULTIPLICITY = 1

# Number of CPU cores to use
CPU = 6

# Amount of memory to use (formats: 5GB, 500MB, 1000MB)
MEMORY = 5GB

# Quantum chemical method (optional, defaults to B3LYP)
METHOD = B3LYP

# Basis set specification (optional)
# If not specified, defaults are: 6-31G* (non-metals), LANL2DZ (metals)
# BASIS_SET = 6-31G*, LANL2DZ

"""
    
    try:
        with open(filename, 'w') as f:
            f.write(sample_content)
        return True
    except Exception as e:
        print(f"[ERROR] Failed to create sample config: {e}")
        return False

def main():
    
    # Check if config file exists, if not create a sample one
    config_file = 'psi4.config'
    
    if not os.path.exists(config_file):
        create_sample_config(config_file)
        sys.exit(0)
    
    # Read configuration
    config = read_config_file(config_file)
    if config is None:
        print("Failed to read configuration file. Exiting.")
        sys.exit(1)
    
    # Read XYZ file
    try:
        xyz_coords = read_xyz_file(config['XYZ'])
    except Exception as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    
    # Run the complete workflow
    energy, optimized_geom, hessian, charges = optimize_metal_complex_with_resp(
        xyz_coordinates=xyz_coords,
        charge=config['CHARGE'],
        multiplicity=config['MULTIPLICITY'],
        memory=config['MEMORY'],
        ncpus=config['CPU'],
        method=config['METHOD'],
        basis_spec=config['BASIS_SET']
    )
    
    if energy is not None:
        if hessian is not None:
            print("  • hessian.txt - Second derivative matrix")
        if charges is not None:
            print("  • resp_charges.dat - RESP atomic charges")
        
        # Create a summary file
        try:
            with open('calculation_summary.txt', 'w') as f:
                f.write("Psi4 Calculation Summary\n")
                f.write("=" * 30 + "\n\n")
                f.write(f"Input file: {config['XYZ']}\n")
                f.write(f"Method: {config['METHOD']}\n")
                f.write(f"Charge: {config['CHARGE']}\n")
                f.write(f"Multiplicity: {config['MULTIPLICITY']}\n")
                if config['BASIS_SET']:
                    f.write(f"Basis set: {config['BASIS_SET']}\n")
                else:
                    f.write(f"Basis set: Default (6-31G*, LANL2DZ)\n")
                f.write(f"Final Energy: {energy:.8f} Hartree\n")
                f.write(f"Memory used: {config['MEMORY']}\n")
                f.write(f"CPUs used: {config['CPU']}\n")
                f.write("\nGenerated Files:\n")
                f.write("• optimized.xyz - Optimized geometry\n")
                if hessian is not None:
                    f.write("• hessian.txt - Hessian matrix\n")
                if charges is not None:
                    f.write("• resp_charges.dat - RESP charges\n")
                f.write("• calculation_summary.txt - This summary\n")
            
        except Exception as e:
            print(f"Warning: Could not create summary file: {e}")
            
    else:
        print("Please check the error messages above and your input files.")

def example_usage():
    main()

if __name__ == "__main__":
    # Run with configuration file
    main()
