#!/bin/bash
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
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


from Bio import PDB
import numpy as np
import re
import shutil
import sys
import argparse
import periodictable as pt

# Import existing functions from second file
def get_atomic_number(element):
    try:
        return getattr(pt, element.lower().capitalize()).number
    except AttributeError:
        base_element = ''.join(c for c in element if not c.isdigit())
        return getattr(pt, base_element.lower().capitalize()).number

# extract the charges from mol2 
def extract_charges_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        charges = []
        is_atom_section = False
        
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                break
            if is_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    charge = float(parts[8])
                    charges.append(charge)
                    
        return charges
    except Exception as e:
        print(f"Error reading charges from {file_path}: {str(e)}")
        return []

def get_perpendicular_vector(v):
    if abs(v[0]) < abs(v[1]):
        perp = np.array([0, -v[2], v[1]])
    else:
        perp = np.array([-v[2], 0, v[0]])
    return perp / np.linalg.norm(perp)

def rotate_vector(v, axis, angle):
    cos_ang = np.cos(angle)
    sin_ang = np.sin(angle)
    return (v * cos_ang + 
            np.cross(axis, v) * sin_ang + 
            axis * np.dot(axis, v) * (1 - cos_ang))

def get_local_environment(residue, target_atom_name, max_bond_distance=1.6):
    target_atom = None
    bonded_atoms = []
    
    for atom in residue:
        if atom.name == target_atom_name:
            target_atom = atom
            break
            
    if target_atom is None:
        return None, []
        
    for atom in residue:
        if atom.name != target_atom_name:
            distance = np.linalg.norm(atom.coord - target_atom.coord)
            if distance <= max_bond_distance:
                bonded_atoms.append((atom, distance))
    
    bonded_atoms.sort(key=lambda x: x[1])
    return target_atom, [atom for atom, _ in bonded_atoms]

def calculate_methyl_position_from_environment(target_atom, bonded_atoms):
    if not bonded_atoms:
        return None, None
        
    target_coord = target_atom.coord
    n_bonds = len(bonded_atoms)
    
    ref_coords = [atom.coord for atom in bonded_atoms]
    
    if target_atom.element == "N":
        if n_bonds == 1:
            ref_vector = ref_coords[0] - target_coord
            ref_vector = ref_vector / np.linalg.norm(ref_vector)
            methyl_direction = -ref_vector
            
        elif n_bonds == 2:
            v1 = ref_coords[0] - target_coord
            v2 = ref_coords[1] - target_coord
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            
            avg_direction = (v1 + v2) / 2
            avg_direction = avg_direction / np.linalg.norm(avg_direction)
            
            methyl_direction = -avg_direction
            
        elif n_bonds == 3:
            v1 = ref_coords[0] - target_coord
            v2 = ref_coords[1] - target_coord
            v3 = ref_coords[2] - target_coord
            v1 = v1 / np.linalg.norm(v1)
            v2 = v2 / np.linalg.norm(v2)
            v3 = v3 / np.linalg.norm(v3)
            
            centroid = (v1 + v2 + v3) / 3
            if np.linalg.norm(centroid) < 1e-6:
                plane_normal = np.cross(v2 - v1, v3 - v1)
                methyl_direction = plane_normal / np.linalg.norm(plane_normal)
            else:
                methyl_direction = -centroid / np.linalg.norm(centroid)
        else:
            return None, None
            
        bond_length = 1.47
        c_pos = target_coord + methyl_direction * bond_length
        h_positions = calculate_methyl_hydrogens(c_pos, target_coord, methyl_direction)
        
    elif target_atom.element == "C":
        avg_direction = np.zeros(3)
        for ref_coord in ref_coords:
            direction = target_coord - ref_coord
            avg_direction += direction / np.linalg.norm(direction)
        avg_direction = avg_direction / n_bonds
        avg_direction = avg_direction / np.linalg.norm(avg_direction)
        
        bond_length = 1.5
        if n_bonds == 1:
            methyl_direction = avg_direction
        elif n_bonds == 2:
            bisector = avg_direction
            plane_normal = np.cross(ref_coords[1] - ref_coords[0], bisector)
            plane_normal = plane_normal / np.linalg.norm(plane_normal)
            methyl_direction = bisector
        elif n_bonds == 3:
            methyl_direction = -avg_direction
        else:
            return None, None
            
        c_pos = target_coord + methyl_direction * bond_length
        h_positions = calculate_methyl_hydrogens(c_pos, target_coord, methyl_direction)
    else:
        return None, None
    
    return c_pos, h_positions

def calculate_methyl_hydrogens(c_pos, target_coord, methyl_direction, h_bond_length=1.09, h_bond_angle=np.radians(109.5)):
    perp1 = get_perpendicular_vector(methyl_direction)
    perp2 = np.cross(methyl_direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    h_positions = []
    rot_angle = np.radians(120)
    
    for i in range(3):
        rotation = i * rot_angle
        h_direction = (np.cos(h_bond_angle) * -methyl_direction + 
                      np.sin(h_bond_angle) * (np.cos(rotation) * perp1 + 
                                            np.sin(rotation) * perp2))
        h_pos = c_pos + h_direction * h_bond_length
        h_positions.append(h_pos)
    
    return h_positions

def calculate_terminal_groups_position(target_atom, bonded_atoms):
    if not bonded_atoms:
        return None
        
    target_coord = target_atom.coord
    n_bonds = len(bonded_atoms)
    
    ref_coords = [atom.coord for atom in bonded_atoms]
    
    if target_atom.element == "N":  # N-terminal: Add acetyl (CO-CH3)
        if n_bonds >= 1:  # As long as N has at least one bond
            ref_vector = ref_coords[0] - target_coord
            ref_vector = ref_vector / np.linalg.norm(ref_vector)
            acetyl_direction = -ref_vector
            
            # Position the C=O carbon (1.47Å from N)
            co_bond_length = 1.47
            co_carbon_pos = target_coord + acetyl_direction * co_bond_length
            
            # Position the O (1.23Å from C=O carbon)
            co_vector = acetyl_direction
            o_pos = co_carbon_pos + co_vector * 1.23
            
            # Position the methyl carbon (1.5Å from C=O carbon)
            ch3_vector = get_perpendicular_vector(co_vector)
            ch3_angle = np.radians(120)  # Tetrahedral angle
            ch3_direction = (np.cos(ch3_angle) * co_vector + 
                           np.sin(ch3_angle) * ch3_vector)
            ch3_pos = co_carbon_pos + ch3_direction * 1.5
            
            # Calculate methyl hydrogens
            h_positions = calculate_methyl_hydrogens(ch3_pos, co_carbon_pos, ch3_direction)
            
            return {
                'type': 'acetyl',
                'positions': {
                    'co_carbon': co_carbon_pos,
                    'oxygen': o_pos,
                    'methyl_carbon': ch3_pos,
                    'hydrogens': h_positions
                }
            }
            
    elif target_atom.element == "C":  # C-terminal: Add NH2
        if n_bonds >= 1:
            avg_direction = np.zeros(3)
            for ref_coord in ref_coords:
                direction = target_coord - ref_coord
                avg_direction += direction / np.linalg.norm(direction)
            avg_direction = avg_direction / n_bonds
            avg_direction = avg_direction / np.linalg.norm(avg_direction)
            
            # Position NH2 group
            n_bond_length = 1.32  # C-N bond length
            n_pos = target_coord + avg_direction * n_bond_length
            
            # Calculate H positions for NH2
            h_bond_length = 1.01  # N-H bond length
            h_angle = np.radians(120)  # H-N-H angle
            
            perp = get_perpendicular_vector(avg_direction)
            h1_direction = (np.cos(h_angle/2) * avg_direction + 
                          np.sin(h_angle/2) * perp)
            h2_direction = (np.cos(h_angle/2) * avg_direction - 
                          np.sin(h_angle/2) * perp)
            
            h1_pos = n_pos + h1_direction * h_bond_length
            h2_pos = n_pos + h2_direction * h_bond_length
            
            return {
                'type': 'amino',
                'positions': {
                    'nitrogen': n_pos,
                    'hydrogen1': h1_pos,
                    'hydrogen2': h2_pos
                }
            }
            
    return None

def get_amide_plane_vectors(target_coord, bonded_atoms):
    # Find alpha carbon and any H atoms
    alpha_c = None
    h_atom = None
    other_heavy = []
    
    for atom in bonded_atoms:
        if atom.name == "CA":
            alpha_c = atom
        elif atom.name.startswith("H"):
            h_atom = atom
        else:
            other_heavy.append(atom)
    
    if alpha_c is None and not other_heavy:
        return None, None
        
    # Use alpha carbon if available, otherwise use first heavy atom
    ref_atom = alpha_c if alpha_c is not None else other_heavy[0]
    ref_vector = target_coord - ref_atom.coord
    ref_vector = ref_vector / np.linalg.norm(ref_vector)
    
    # Get perpendicular vector avoiding H if present
    if h_atom:
        h_vector = h_atom.coord - target_coord
        h_vector = h_vector / np.linalg.norm(h_vector)
        plane_normal = np.cross(ref_vector, h_vector)
    else:
        # Use any valid perpendicular vector
        plane_normal = get_perpendicular_vector(ref_vector)
    
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    
    return plane_normal, ref_vector
def calculate_acetyl_position_from_environment(target_atom, bonded_atoms):
    if not bonded_atoms or len(bonded_atoms) < 1:
        return None, None, None
        
    target_coord = target_atom.coord
    
    # Get amide plane geometry
    plane_normal, ca_n_vector = get_amide_plane_vectors(target_coord, bonded_atoms)
    if plane_normal is None:
        return None, None, None
    
    # Calculate C-N bond direction in the amide plane
    # Should be roughly opposite to CA-N but slightly rotated
    c_n_angle = np.radians(123)  # Typical C-N-CA angle in amides
    c_n_vector = rotate_vector(-ca_n_vector, plane_normal, c_n_angle)
    
    # Place carbonyl carbon
    c_n_bond = 1.335  # Amide C-N bond length
    c_pos = target_coord + c_n_vector * c_n_bond
    
    # Calculate carbonyl oxygen position
    # O=C-N angle should be ~121Â° and in the amide plane
    o_c_n_angle = np.radians(121)
    o_direction = rotate_vector(-c_n_vector, plane_normal, o_c_n_angle)
    c_o_bond = 1.229  # C=O bond length
    o_pos = c_pos + o_direction * c_o_bond
    
    # Calculate methyl carbon position
    # Should be roughly tetrahedral relative to C=O and C-N bonds
    c_c_bond = 1.508  # C-C single bond length
    ch3_c_n_angle = np.radians(85)  # Typical CH3-C-N angle
    
    # Rotate in opposite direction from oxygen to maintain proper geometry
    ch3_direction = rotate_vector(-c_n_vector, plane_normal, -ch3_c_n_angle)
    ch3_pos = c_pos + ch3_direction * c_c_bond
    
    return c_pos, o_pos, ch3_pos

def calculate_acetyl_group(n_atom, bonded_atoms):
    if not bonded_atoms or len(bonded_atoms) == 0:
        return None, None, None, None
        
    # Calculate core acetyl positions
    c_pos, o_pos, ch3_pos = calculate_acetyl_position_from_environment(n_atom, bonded_atoms)
    
    if c_pos is None:
        return None, None, None, None
    
    # Calculate methyl hydrogen positions ensuring proper tetrahedral geometry
    c_ch3_vector = ch3_pos - c_pos
    c_ch3_vector = c_ch3_vector / np.linalg.norm(c_ch3_vector)
    
    # Calculate reference vector for hydrogen placement
    # Use C=O bond as reference to ensure proper staggering
    co_vector = o_pos - c_pos
    co_vector = co_vector / np.linalg.norm(co_vector)
    
    # Calculate hydrogens with 109.5Â° tetrahedral angles
    # and proper staggered conformation relative to C=O
    perp1 = np.cross(c_ch3_vector, co_vector)
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(c_ch3_vector, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    h_bond_length = 1.090  # C-H bond length
    h_positions = []
    
    for i in range(3):
        angle = i * np.radians(120)
        h_direction = (np.cos(np.radians(109.5)) * -c_ch3_vector +
                      np.sin(np.radians(109.5)) * 
                      (np.cos(angle) * perp1 + np.sin(angle) * perp2))
        h_pos = ch3_pos + h_direction * h_bond_length
        h_positions.append(h_pos)
    
    return c_pos, o_pos, ch3_pos, h_positions

def calculate_amino_group(c_atom, bonded_atoms):
    if not bonded_atoms or len(bonded_atoms) == 0:
        return None, None
        
    c_coord = c_atom.coord
    n_bonds = len(bonded_atoms)
    
    # Calculate direction for NH2 attachment
    avg_direction = np.zeros(3)
    for atom in bonded_atoms:
        direction = atom.coord - c_coord
        direction = direction / np.linalg.norm(direction)
        avg_direction += direction
    
    avg_direction = -avg_direction / n_bonds
    avg_direction = avg_direction / np.linalg.norm(avg_direction)
    
    # Place nitrogen
    c_n_bond = 1.35  # C-N bond length
    n_pos = c_coord + avg_direction * c_n_bond
    
    # Place hydrogens
    n_h_bond = 1.01  # N-H bond length
    perp = get_perpendicular_vector(avg_direction)
    h_positions = []
    
    for angle in [-60, 60]:
        h_direction = rotate_vector(avg_direction, perp, np.radians(angle))
        h_pos = n_pos + h_direction * n_h_bond
        h_positions.append(h_pos)
    
    return n_pos, h_positions

#Identify N-terminal and C-terminal residues for each chain.
def identify_terminal_residues(structure, standard_residues):
    
    terminal_info = {}
    
    # Base residue names (without N/C prefix)
    base_residues = {res for res in standard_residues 
                     if not res.startswith(('N', 'C')) and res not in {'NHE', 'NME', 'ACE'}}
    
    for model in structure:
        for chain in model:
            residues = [res for res in chain if res.get_resname() in base_residues]
            
            if len(residues) == 0:
                continue
            
            # Mark first standard residue as N-terminal
            first_res = residues[0]
            terminal_info[(chain.id, first_res.get_id())] = 'NTERM'
            
            # Mark last standard residue as C-terminal
            last_res = residues[-1]
            terminal_info[(chain.id, last_res.get_id())] = 'CTERM'
    
    return terminal_info

def analyze_and_extract_metal_site(input_pdb, mol2_file, metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5'], distance_cutoff=2.6):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    if isinstance(metals, list):
        metals = set(metals)
    elif not isinstance(metals, set):
        metals = {'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU',
                'NA', 'K', 'CA', 'LI', 'RB', 'CS', 'MG', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI',
                'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}

    standard_residues = {
            # Base residues
            "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
            "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "NHE", "NME",
            "ACE",
            # N-terminal variants
            "NALA", "NARG", "NASH", "NASN", "NASP", "NCYM", "NCYS", "NCYX", "NGLH", "NGLN",
            "NGLU", "NGLY", "NHID", "NHIE", "NHIP", "NHYP", "NILE", "NLEU", "NLYN", "NLYS",
            "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL",
            # C-terminal variants
            "CALA", "CARG", "CASH", "CASN", "CASP", "CCYM", "CCYS", "CCYX", "CGLH", "CGLN",
            "CGLU", "CGLY", "CHID", "CHIE", "CHIP", "CHYP", "CILE", "CLEU", "CLYN", "CLYS",
            "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL"
    }

    # Identify terminals
    terminal_info = identify_terminal_residues(structure, standard_residues)

    mol2_charges = extract_charges_from_mol2(mol2_file)
    if not mol2_charges:
        raise Exception("Failed to extract charges from MOL2 file")

    metal_coordination = {}
    residues_to_extract = set()
    coordinated_standard_residues = set()
    original_order = []
    atoms_data = []
    coordinated_residue_counts = {}

    # Analysis phase
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_key = (chain.id, residue.get_id())
                original_order.append(residue_key)

                for atom in residue:
                    if atom.element in metals:
                        metal_coord = atom.coord
                        key = f"{atom.element}_{chain.id}_{residue.get_id()[1]}"
                        metal_coordination[key] = {
                            'metal_position': metal_coord,
                            'metal_element': atom.element,
                            'coordinating_residues': [],
                            'coordination_number': 0
                        }

                        for chain2 in model:
                            for residue2 in chain2:
                                is_coordinating = False
                                coordinating_atoms = []

                                for atom2 in residue2:
                                    distance = np.linalg.norm(metal_coord - atom2.coord)
                                    if distance <= distance_cutoff:
                                        is_coordinating = True
                                        coord_info = {
                                            'chain': chain2.id,
                                            'residue_number': residue2.get_id()[1],
                                            'residue_name': residue2.get_resname(),
                                            'atom_name': atom2.get_name().strip(),
                                            'distance': round(distance, 2),
                                            'atom_coord': atom2.coord,
                                            'element': atom2.element
                                        }
                                        coordinating_atoms.append(coord_info)
                                        metal_coordination[key]['coordination_number'] += 1

                                if is_coordinating:
                                    metal_coordination[key]['coordinating_residues'].extend(coordinating_atoms)
                                    residues_to_extract.add((chain2.id, residue2.get_id()))

                                    if residue2.get_resname() in standard_residues:
                                        coordinated_standard_residues.add((chain2.id, residue2.get_id()))
                                        res_name = residue2.get_resname()
                                        coordinated_residue_counts[res_name] = coordinated_residue_counts.get(res_name, 0) + 1

    # Find standard residues linked to non-standard coordinating residues
    heavy_atoms = [atom for model in structure for chain in model for residue in chain for atom in residue if atom.element != 'H']
    ns = PDB.NeighborSearch(heavy_atoms)
    bond_cutoff = 1.9

    nonstandard_coordinating = [(chain_id, res_id) for chain_id, res_id in residues_to_extract
                                if structure[0][chain_id][res_id].resname not in standard_residues]

    for chain_id, res_id in nonstandard_coordinating:
        nonstandard_res = structure[0][chain_id][res_id]
        nonstandard_heavy_atoms = [atom for atom in nonstandard_res if atom.element != 'H']

        for atom in nonstandard_heavy_atoms:
            nearby_atoms = ns.search(atom.coord, bond_cutoff, level='A')
            for nearby_atom in nearby_atoms:
                nearby_res = nearby_atom.get_parent()
                nearby_chain = nearby_res.get_parent()

                if (nearby_chain.id != chain_id or nearby_res.get_id() != res_id) and nearby_res.resname in standard_residues:
                    nearby_key = (nearby_chain.id, nearby_res.get_id())
                    if nearby_key not in residues_to_extract:
                        residues_to_extract.add(nearby_key)
                        if nearby_res.resname in standard_residues:
                            coordinated_standard_residues.add(nearby_key)
                            res_name = nearby_res.resname
                            coordinated_residue_counts[res_name] = coordinated_residue_counts.get(res_name, 0) + 1

    # Improved cross-residue bond detection
    extracted_atoms = []
    extracted_residue_objects = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                res_key = (chain.id, residue.get_id())
                if res_key in residues_to_extract:
                    extracted_residue_objects[res_key] = residue
                    for atom in residue:
                        if atom.element != 'H':
                            extracted_atoms.append(atom)

    ns_extracted = PDB.NeighborSearch(extracted_atoms)

    cross_residue_bonds = {}
    for res_key in residues_to_extract:
        cross_residue_bonds[res_key] = {}

    for atom in extracted_atoms:
        res_key = (atom.get_parent().get_parent().id, atom.get_parent().get_id())
        atom_name = atom.name

        if atom_name not in cross_residue_bonds[res_key]:
            cross_residue_bonds[res_key][atom_name] = []

        nearby_atoms = ns_extracted.search(atom.coord, bond_cutoff, level='A')
        for nearby_atom in nearby_atoms:
            if nearby_atom == atom:
                continue

            nearby_res_key = (nearby_atom.get_parent().get_parent().id,
                             nearby_atom.get_parent().get_id())

            if nearby_res_key != res_key:
                cross_residue_bonds[res_key][atom_name].append(nearby_atom)

    # Lists to store all atoms (original and capping)
    all_atoms = []
    terminal_modifications = []

    # Process residues for extraction
    for chain_id, res_id in original_order:
        if (chain_id, res_id) in residues_to_extract:
            residue = structure[0][chain_id][res_id]
            is_standard = residue.resname in standard_residues

            # Add original atoms
            for atom in residue:
                element = atom.element if atom.element != " " else atom.name[0]
                all_atoms.append({
                    'element': element,
                    'coord': atom.coord,
                    'name': atom.name,
                    'resname': residue.resname,
                    'resnum': res_id[1],
                    'is_capping': False
                })
                atoms_data.append({
                    'atomic_number': get_atomic_number(element),
                    'is_standard': 0
                })

            # Add terminal groups for standard residues that coordinate with metals
            if (chain_id, res_id) in coordinated_standard_residues:
                res_key = (chain_id, res_id)

                # Check if N atom needs acetylation
                n_atom, n_bonded = get_local_environment(residue, "N", max_bond_distance=1.6)
                if n_atom is not None:
                    needs_acetyl = True

                    if "N" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["N"]:
                        needs_acetyl = False

                    if len(n_bonded) >= 3:
                        needs_acetyl = False

                    if needs_acetyl:
                        c_pos, o_pos, ch3_pos, h_positions = calculate_acetyl_group(n_atom, n_bonded)
                        if c_pos is not None:
                            # Add C=O carbon
                            terminal_modifications.append({
                                'element': 'C', 'coord': c_pos, 'name': 'CAC',
                                'resname': residue.resname, 'resnum': res_id[1],
                                'is_capping': True
                            })
                            # Add O
                            terminal_modifications.append({
                                'element': 'O', 'coord': o_pos, 'name': 'OAC',
                                'resname': residue.resname, 'resnum': res_id[1],
                                'is_capping': True
                            })
                            # Add methyl carbon
                            terminal_modifications.append({
                                'element': 'C', 'coord': ch3_pos, 'name': 'CME',
                                'resname': residue.resname, 'resnum': res_id[1],
                                'is_capping': True
                            })
                            # Add methyl hydrogens
                            for i, h_pos in enumerate(h_positions):
                                terminal_modifications.append({
                                    'element': 'H', 'coord': h_pos, 'name': f'HM{i+1}',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                })

                # Check if C atom needs amination
                c_atom, c_bonded = get_local_environment(residue, "C", max_bond_distance=1.6)
                if c_atom is not None:
                    needs_amination = True

                    if "C" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["C"]:
                        needs_amination = False

                    if len(c_bonded) >= 3:
                        needs_amination = False

                    if needs_amination:
                        n_pos, h_positions = calculate_amino_group(c_atom, c_bonded)
                        if n_pos is not None:
                            # Add N
                            terminal_modifications.append({
                                'element': 'N', 'coord': n_pos, 'name': 'NT',
                                'resname': residue.resname, 'resnum': res_id[1],
                                'is_capping': True
                            })
                            # Add H atoms
                            for i, h_pos in enumerate(h_positions):
                                terminal_modifications.append({
                                    'element': 'H', 'coord': h_pos, 'name': f'HT{i+1}',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                })

    # Add terminal modifications to all_atoms
    all_atoms.extend(terminal_modifications)

    # Add terminal modifications to atoms_data for charge calculations
    for mod in terminal_modifications:
        atoms_data.append({
            'atomic_number': get_atomic_number(mod['element']),
            'is_standard': 1
        })

    # Write reference_structure.xyz
    with open('reference_structure.xyz', 'w') as f:
        f.write(f"{len(all_atoms)}\n\n")
        for atom in all_atoms:
            if atom['is_capping']:
                f.write(f"{atom['element']:2} {atom['coord'][0]:10.6f} {atom['coord'][1]:10.6f} {atom['coord'][2]:10.6f} {atom['name']:4} {atom['resname']:3} {atom['resnum']:4}  CAPPING\n")
            else:
                f.write(f"{atom['element']:2} {atom['coord'][0]:10.6f} {atom['coord'][1]:10.6f} {atom['coord'][2]:10.6f} {atom['name']:4} {atom['resname']:3} {atom['resnum']:4}\n")

    # Write processed charges
    with open('processed_charges.dat', 'w') as charge_file:
        if len(atoms_data) != len(mol2_charges):
            raise Exception(f"Number of atoms mismatch: Structure has {len(atoms_data)}, MOL2 has {len(mol2_charges)}")

        for i, atom_data in enumerate(atoms_data):
            charge = 0.000000 if atom_data['is_standard'] == 1 else mol2_charges[i]
            charge_file.write(f"{charge:.6f}\n")

    # Write capping_link_atoms.dat
    # Read reference_structure.xyz as the single source of truth for atom ordering.
    # Skip the two header lines, count every atom line (1-based), and record the
    # index of any CAPPING atom whose name is CAC (acetyl carbonyl C) or NT (amino N).
    with open('reference_structure.xyz', 'r') as xyz, \
         open('capping_link_atoms.dat', 'w') as out:

        atom_index = 0
        for line in xyz:
            parts = line.split()
            # Skip header lines: the count line (single integer) and the blank comment line
            if len(parts) < 6:
                continue
            if parts[0].lstrip('-').replace('.', '', 1).isdigit():
                continue

            atom_index += 1
            atom_name  = parts[4]           # column 5: CAC, OAC, CME, NT, HT1 ...
            is_capping = 'CAPPING' in parts  # last column present only on capping lines

            if is_capping and atom_name == 'CAC':
                out.write(f"{atom_index} C\n")
            elif is_capping and atom_name == 'NT':
                out.write(f"{atom_index} N\n")

    # Generate summary of coordinated residues
    coordinated_residue_list = [(res_name, count) for res_name, count in coordinated_residue_counts.items()]

    # Write coordination analysis
    with open('coordination_analysis.txt', 'w') as f:
        f.write("Metal Coordination Analysis\n")
        f.write("==========================\n\n")

        for metal_key, info in metal_coordination.items():
            f.write(f"Metal: {info['metal_element']} (Coordination number: {info['coordination_number']})\n")
            f.write("Coordinating residues:\n")

            for res in info['coordinating_residues']:
                f.write(f"  {res['residue_name']} {res['residue_number']} "
                       f"(Chain {res['chain']}) - {res['atom_name']} "
                       f"at {res['distance']} Å\n")
            f.write("\n")

        f.write("\nCAPPING Group Additions:\n")
        f.write("=====================\n")
        acetyl_count = sum(1 for atom in terminal_modifications if atom['name'] in ['CAC', 'CME'])
        amino_count  = sum(1 for atom in terminal_modifications if atom['name'] == 'NT')
        f.write(f"Total acetyl groups added: {acetyl_count//2}\n")
        f.write(f"Total amino groups added: {amino_count}\n")
        f.write(f"Total additional atoms added: {len(terminal_modifications)}\n\n")

    with open('terminal_info.dat', 'w') as f:
        for (chain_id, res_id), terminal_type in terminal_info.items():
            f.write(f"{chain_id} {res_id[1]} {terminal_type}\n")

    return metal_coordination, coordinated_residue_list, terminal_info

#Provides default charge for N,C,O and CA atoms from amber force field.
def get_default_reference_charges():
    return {
    'ALA': {'N': -0.4157, 'CA': 0.0337, 'C': 0.5973, 'O': -0.5679},
    'ARG': {'N': -0.3479, 'CA': -0.2637, 'C': 0.7341, 'O': -0.5894},
    'ASH': {'N': -0.4157, 'CA': 0.0341, 'C': 0.5973, 'O': -0.5679},
    'ASN': {'N': -0.4157, 'CA': 0.0143, 'C': 0.5973, 'O': -0.5679},
    'ASP': {'N': -0.5163, 'CA': 0.0381, 'C': 0.5366, 'O': -0.5819},
    'CYM': {'N': -0.4157, 'CA': -0.0351, 'C': 0.5973, 'O': -0.5679},
    'CYS': {'N': -0.4157, 'CA': 0.0213, 'C': 0.5973, 'O': -0.5679},
    'CYX': {'N': -0.4157, 'CA': 0.0429, 'C': 0.5973, 'O': -0.5679},
    'GLH': {'N': -0.4157, 'CA': 0.0145, 'C': 0.5973, 'O': -0.5679},
    'GLN': {'N': -0.4157, 'CA': -0.0031, 'C': 0.5973, 'O': -0.5679},
    'GLU': {'N': -0.5163, 'CA': 0.0397, 'C': 0.5366, 'O': -0.5819},
    'GLY': {'N': -0.4157, 'CA': -0.0252, 'C': 0.5973, 'O': -0.5679},
    'HID': {'N': -0.4157, 'CA': 0.0188, 'C': 0.5973, 'O': -0.5679},
    'HIE': {'N': -0.4157, 'CA': -0.0581, 'C': 0.5973, 'O': -0.5679},
    'HIP': {'N': -0.3479, 'CA': -0.1354, 'C': 0.7341, 'O': -0.5894},
    'HYP': {'N': -0.2548, 'CA': 0.0047, 'C': 0.5896, 'O': -0.5748},
    'ILE': {'N': -0.4157, 'CA': -0.0597, 'C': 0.5973, 'O': -0.5679},
    'LEU': {'N': -0.4157, 'CA': -0.0518, 'C': 0.5973, 'O': -0.5679},
    'LYN': {'N': -0.4157, 'CA': -0.07206, 'C': 0.5973, 'O': -0.5679},
    'LYS': {'N': -0.3479, 'CA': -0.24, 'C': 0.7341, 'O': -0.5894},
    'MET': {'N': -0.4157, 'CA': -0.0237, 'C': 0.5973, 'O': -0.5679},
    'PHE': {'N': -0.4157, 'CA': -0.0024, 'C': 0.5973, 'O': -0.5679},
    'PRO': {'N': -0.2548, 'CA': -0.0266, 'C': 0.5896, 'O': -0.5748},
    'SER': {'N': -0.4157, 'CA': -0.0249, 'C': 0.5973, 'O': -0.5679},
    'THR': {'N': -0.4157, 'CA': -0.0389, 'C': 0.5973, 'O': -0.5679},
    'TRP': {'N': -0.4157, 'CA': -0.0275, 'C': 0.5973, 'O': -0.5679},
    'TYR': {'N': -0.4157, 'CA': -0.0014, 'C': 0.5973, 'O': -0.5679},
    'VAL': {'N': -0.4157, 'CA': -0.0875, 'C': 0.5973, 'O': -0.5679},
    'CALA': {'N': -0.3821, 'CA': -0.1747, 'C': 0.7731, 'O': -0.8055},
    'CARG': {'N': -0.3481, 'CA': -0.3068, 'C': 0.8557, 'O': -0.8266},
    'CASN': {'N': -0.3821, 'CA': -0.208, 'C': 0.805, 'O': -0.8147},
    'CASP': {'N': -0.5192, 'CA': -0.1817, 'C': 0.7256, 'O': -0.7887},
    'CCYS': {'N': -0.3821, 'CA': -0.1635, 'C': 0.7497, 'O': -0.7981},
    'CCYX': {'N': -0.3821, 'CA': -0.1318, 'C': 0.7618, 'O': -0.8041},
    'CGLN': {'N': -0.3821, 'CA': -0.2248, 'C': 0.7775, 'O': -0.8042},
    'CGLU': {'N': -0.5192, 'CA': -0.2059, 'C': 0.742, 'O': -0.793},
    'CGLY': {'N': -0.3821, 'CA': -0.2493, 'C': 0.7231, 'O': -0.7855},
    'CHID': {'N': -0.3821, 'CA': -0.1739, 'C': 0.7615, 'O': -0.8016},
    'CHIE': {'N': -0.3821, 'CA': -0.2699, 'C': 0.7916, 'O': -0.8065},
    'CHIP': {'N': -0.3481, 'CA': -0.1445, 'C': 0.8032, 'O': -0.8177},
    'CHYP': {'N': -0.2802, 'CA': -0.0993, 'C': 0.6631, 'O': -0.7697},
    'CILE': {'N': -0.3821, 'CA': -0.31, 'C': 0.8343, 'O': -0.819},
    'CLEU': {'N': -0.3821, 'CA': -0.2847, 'C': 0.8326, 'O': -0.8199},
    'CLYS': {'N': -0.3481, 'CA': -0.2903, 'C': 0.8488, 'O': -0.8252},
    'CMET': {'N': -0.3821, 'CA': -0.2597, 'C': 0.8013, 'O': -0.8105},
    'CPHE': {'N': -0.3821, 'CA': -0.1825, 'C': 0.766, 'O': -0.8026},
    'CPRO': {'N': -0.2802, 'CA': -0.1336, 'C': 0.6631, 'O': -0.7697},
    'CSER': {'N': -0.3821, 'CA': -0.2722, 'C': 0.8113, 'O': -0.8132},
    'CTHR': {'N': -0.3821, 'CA': -0.242, 'C': 0.781, 'O': -0.8044},
    'CTRP': {'N': -0.3821, 'CA': -0.2084, 'C': 0.7658, 'O': -0.8011},
    'CTYR': {'N': -0.3821, 'CA': -0.2015, 'C': 0.7817, 'O': -0.807},
    'CVAL': {'N': -0.3821, 'CA': -0.3438, 'C': 0.835, 'O': -0.8173},
    'NHE': {'N': -0.463},
    'NME': {'N': -0.4157, 'C': -0.149},
    'ACE': {'C': 0.5972, 'O': -0.5679},
    'NALA': {'N': 0.1414, 'CA': 0.0962, 'C': 0.6163, 'O': -0.5722},
    'NARG': {'N': 0.1305, 'CA': -0.0223, 'C': 0.7214, 'O': -0.6013},
    'NASN': {'N': 0.1801, 'CA': 0.0368, 'C': 0.6163, 'O': -0.5722},
    'NASP': {'N': 0.0782, 'CA': 0.0292, 'C': 0.5621, 'O': -0.5889},
    'NCYS': {'N': 0.1325, 'CA': 0.0927, 'C': 0.6123, 'O': -0.5713},
    'NCYX': {'N': 0.2069, 'CA': 0.1055, 'C': 0.6123, 'O': -0.5713},
    'NGLN': {'N': 0.1493, 'CA': 0.0536, 'C': 0.6123, 'O': -0.5713},
    'NGLU': {'N': 0.0017, 'CA': 0.0588, 'C': 0.5621, 'O': -0.5889},
    'NGLY': {'N': 0.2943, 'CA': -0.01, 'C': 0.6163, 'O': -0.5722},
    'NHID': {'N': 0.1542, 'CA': 0.0964, 'C': 0.6123, 'O': -0.5713},
    'NHIE': {'N': 0.1472, 'CA': 0.0236, 'C': 0.6123, 'O': -0.5713},
    'NHIP': {'N': 0.256, 'CA': 0.0581, 'C': 0.7214, 'O': -0.6013},
    'NILE': {'N': 0.0311, 'CA': 0.0257, 'C': 0.6123, 'O': -0.5713},
    'NLEU': {'N': 0.101, 'CA': 0.0104, 'C': 0.6123, 'O': -0.5713},
    'NLYS': {'N': 0.0966, 'CA': -0.0015, 'C': 0.7214, 'O': -0.6013},
    'NMET': {'N': 0.1592, 'CA': 0.0221, 'C': 0.6123, 'O': -0.5713},
    'NPHE': {'N': 0.1737, 'CA': 0.0733, 'C': 0.6123, 'O': -0.5713},
    'NPRO': {'N': -0.202, 'CA': 0.1, 'C': 0.526, 'O': -0.5},
    'NSER': {'N': 0.1849, 'CA': 0.0567, 'C': 0.6163, 'O': -0.5722},
    'NTHR': {'N': 0.1812, 'CA': 0.0034, 'C': 0.6163, 'O': -0.5722},
    'NTRP': {'N': 0.1913, 'CA': 0.0421, 'C': 0.6123, 'O': -0.5713},
    'NTYR': {'N': 0.194, 'CA': 0.057, 'C': 0.6123, 'O': -0.5713},
    'NVAL': {'N': 0.0577, 'CA': -0.0054, 'C': 0.6163, 'O': -0.5722}}

#Read the fixed charges file 
def read_fixed_charges_file(filename, terminal_info_file='terminal_info.dat'):
    # Read terminal info
    terminal_lookup = {}
    try:
        with open(terminal_info_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    chain_id = parts[0]
                    resnum = int(parts[1])
                    terminal_type = parts[2]
                    terminal_lookup[resnum] = terminal_type
    except FileNotFoundError:
        print(f"Warning: {terminal_info_file} not found. No terminal info available.")

    # Read fixed charges
    fixed_charges_map = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    atom_id = int(parts[0])
                    atom_name_ref = parts[1]
                    atom_name_out = parts[2]
                    residue_name = parts[3]

                    fixed_charges_map[(atom_id, atom_name_out)] = (atom_name_ref, residue_name)

    except FileNotFoundError:
        print(f"Warning: {filename} not found. Proceeding without fixed charges.")
        return {}, terminal_lookup

    return fixed_charges_map, terminal_lookup

#Process charges while using reference charges from get_default_reference_charges for N, C, O and CA atoms.
def process_charges(processed_charges_file, mol2_file, target_charge, coordinated_residues, fixed_charges_file):
    # Get reference charges
    reference_charges = get_default_reference_charges()
    
    # Read the fixed charges file
    fixed_charges_map, terminal_lookup = read_fixed_charges_file(fixed_charges_file)
 
    # Read the original charges
    with open(processed_charges_file, 'r') as f:
        charges = [float(line.strip()) for line in f]
    
    # Make a copy of original charges
    new_charges = charges.copy()
    
    # First handle the zero charges distribution
    # Get non-zero charges and their indices
    zero_indices = set(i for i, charge in enumerate(charges) if abs(charge) <= 1e-6)
    non_zero_indices = [i for i, charge in enumerate(charges) if abs(charge) > 1e-6]
    non_zero_charges = [charges[i] for i in non_zero_indices]
    total_number_atoms = len(non_zero_indices)
    
    # Calculate total charge of non-zero elements
    total_selected_charges = sum(non_zero_charges)
    
    # Calculate correction term for zero charges
    correction_term = (target_charge - total_selected_charges) / total_number_atoms
    
    # Distribute correction among non-zero charges
    for idx in non_zero_indices:
        new_charges[idx] += correction_term
    
    # Now proceed with MOL2 file reading and fixed charges
    is_atom_section = False
    atom_index = 0
    fixed_indices = set()  # Track atoms with reference charges
    mol2_atoms = []  # Store all atoms with their details
    
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                is_atom_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                is_atom_section = False
                continue
            
            if is_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    mol2_atom_id = int(parts[0])
                    atom_type = parts[5]
                    residue_name = parts[7]
                    residue_number = int(parts[6])  # Get residue number from MOL2
                
                    # Store atom information
                    mol2_atoms.append({
                        'index': atom_index,
                        'id': mol2_atom_id,
                        'type': atom_type,
                        'residue': residue_name,
                        'resnum': residue_number

                    })
                    
                    # Check if this atom should have a reference charge
                    key = (mol2_atom_id, atom_type)
                    if key in fixed_charges_map:
                        atom_name_ref, ref_residue = fixed_charges_map[key]
                        
                        # Check if this residue is a terminal
                        lookup_residue = ref_residue
                        if residue_number in terminal_lookup:
                            terminal_type = terminal_lookup[residue_number]
                            # Add N or C prefix for terminal lookup
                            lookup_residue = ('N' if terminal_type == 'NTERM' else 'C') + ref_residue
                        
                        # Look up charge using terminal-aware residue name
                        if lookup_residue in reference_charges and atom_name_ref in reference_charges[lookup_residue]:
                            new_charges[atom_index] = reference_charges[lookup_residue][atom_name_ref]
                            fixed_indices.add(atom_index)

                    atom_index += 1
    
    # Calculate total fixed charge after applying reference charges
    fixed_total = sum(new_charges[idx] for idx in fixed_indices)
    # Calculate remaining charge to distribute
    remaining_charge = target_charge - fixed_total
    # Get indices of non-fixed atoms, excluding both fixed and zero indices
    non_fixed_indices = [i for i in range(len(charges)) 
                        if i not in fixed_indices and i not in zero_indices]
    
    if non_fixed_indices:
        # Calculate original total of non-fixed charges
        original_non_fixed_total = sum(new_charges[i] for i in non_fixed_indices)
        # Calculate scaling factor based on number of non-fixed atoms
        scaling_factor = (remaining_charge - original_non_fixed_total) / len(non_fixed_indices)
        # Apply scaling factor to non-fixed charges
        for idx in non_fixed_indices:
            new_charges[idx] += scaling_factor
    
    # Write output files
    write_charge_files(charges, new_charges, len(non_fixed_indices), 
                      sum(charges), target_charge, remaining_charge,
                      coordinated_residues)
    
    return new_charges
#Write detailed charge information to files
def write_charge_files(original_charges, new_charges, total_number_atoms, 
                      total_selected_charges, target_charge, right_charge,
                      coordinated_residues):
    # Write summary statistics
    with open('charge_statistics.txt', 'w') as f:
        f.write("Charge Distribution Analysis\n")
        f.write("===========================\n\n")
        f.write("Coordinated Standard Residues:\n")
        for residue_name, count in coordinated_residues:
            f.write(f"  {residue_name}: {count} residue(s)\n")
        f.write(f"Target system charge: {target_charge:.6f}\n")
        f.write(f"Charge to distribute: {right_charge:.6f}\n")
        f.write(f"\nDistribution Statistics:\n")
        f.write(f"  Total number of atoms for charge distribution: {total_number_atoms}\n")
        f.write(f"  Sum of original charges: {total_selected_charges:.6f}\n")
        f.write(f"  Correction term per atom: {(right_charge - total_selected_charges) / total_number_atoms:.6f}\n")
    
    # Write new charges
    with open('recalculated_charges.dat', 'w') as f:
        for charge in new_charges:
            f.write(f"{charge:.6f}\n")
    
    # Write original charges
    with open('original_charges.dat', 'w') as f:
        for charge in original_charges:
            f.write(f"{charge:.6f}\n")


#Update charges in a MOL2 file using values from a charge file.    
def update_mol2_file(input_mol2, charge_file='recalculated_charges.dat', output_mol2='updated_NEW_COMPLEX.mol2'):
    try:
        # Read new charges
        new_charges = []
        with open(charge_file, 'r') as f:
            for line in f:
                try:
                    charge = float(line.strip())
                    new_charges.append(charge)
                except ValueError:
                    continue

        # Process MOL2 file
        updated_mol2_lines = []
        current_atom_index = 0
        is_atom_section = False
        
        with open(input_mol2, 'r') as input_file:
            for line in input_file:
                # Check for section headers
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    continue
                elif line.startswith("@<TRIPOS>"):
                    is_atom_section = False
                    updated_mol2_lines.append(line)
                    continue

                # Process atom section
                if is_atom_section and line.strip():
                    try:
                        parts = line.split()
                        if len(parts) >= 9:  # MOL2 format should have at least 9 columns
                            # Format each column with proper spacing
                            atom_id = f"{parts[0]:>7}"                    # Atom ID
                            atom_name = f"{parts[1]:<3}"                  # Atom name
                            x = f"{float(parts[2]):>10.4f}"              # X coordinate
                            y = f"{float(parts[3]):>10.4f}"              # Y coordinate
                            z = f"{float(parts[4]):>10.4f}"              # Z coordinate
                            atom_type = f"{parts[5]:<5}"                 # Atom type
                            subst_id = f"{parts[6]:4>}"                  # Substructure ID
                            subst_name = f"{parts[7]:<4}"                # Substructure name
                            
                            # Add new charge
                            if current_atom_index < len(new_charges):
                                charge = new_charges[current_atom_index]
                                current_atom_index += 1
                            else:
                                charge = float(parts[8])
                            
                            charge_str = f"{charge:>10.6f}"              # Charge
                            
                            # Combine all parts with proper spacing
                            updated_line = f"{atom_id} {atom_name}  {x}  {y}  {z} {atom_type}  {subst_id} {subst_name}   {charge_str}\n"
                            updated_mol2_lines.append(updated_line)
                        else:
                            updated_mol2_lines.append(line)
                    except (ValueError, IndexError):
                        updated_mol2_lines.append(line)
                else:
                    updated_mol2_lines.append(line)

        # Write the updated MOL2 file
        with open(output_mol2, 'w') as output_file:
            output_file.writelines(updated_mol2_lines)

    except FileNotFoundError as e:
        print(f"Error: File not found - {str(e)}")
        raise
    except PermissionError:
        print("Error: Permission denied when accessing files")
        raise
    except Exception as e:
        print(f"Unexpected error occurred: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Analyze metal coordination sites ')
    parser.add_argument('pdb_file', help='Input PDB file to analyze')
    parser.add_argument('mol2_file', help='Input MOL2 file containing charges')
    parser.add_argument('target_charge', type=float, help='Target total charge for the system')

    args = parser.parse_args()

    try:
        metal_sites, coordinated_residues, terminal_info = analyze_and_extract_metal_site(
            args.pdb_file,
            mol2_file=args.mol2_file,
        )
        
        process_charges('processed_charges.dat', args.mol2_file, args.target_charge, coordinated_residues, fixed_charges_file="fixed_charges.dat") 
        # Update MOL2 file with new charges
        update_mol2_file(args.mol2_file, 'recalculated_charges.dat', 'updated_easy_COMPLEX.mol2')
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {str(e)}")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing files: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
