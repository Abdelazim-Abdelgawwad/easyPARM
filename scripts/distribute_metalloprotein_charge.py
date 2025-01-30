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
#                              |  $$$$$$/              Ver. 3.00 - 12 January 2024                                #
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

def analyze_and_extract_metal_site(input_pdb, mol2_file, distance_cutoff=2.6):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)
    
    metals = {'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 
              'NA', 'K', 'CA', 'LI', 'RB', 'CS', 'MG', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 
              'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}

    standard_residues = {
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
        "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    }

    mol2_charges = extract_charges_from_mol2(mol2_file)
    if not mol2_charges:
        raise Exception("Failed to extract charges from MOL2 file")

    metal_coordination = {}
    residues_to_extract = set()
    coordinated_standard_residues = set()
    original_order = []
    atoms_data = []
    coordinated_residue_counts = {}
    
    # Lists to store original and added atoms
    original_atoms = []
    added_groups = []

    # Analysis phase (same as before)
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_key = (chain.id, residue.get_id())
                original_order.append(residue_key)

                for atom in residue:
                    if atom.element in metals:
                        metal_coord = atom.coord
                        key = f"{atom.element}_{chain.id}_{residue.get_id()}"
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

    for chain_id, res_id in original_order:
        if (chain_id, res_id) in residues_to_extract:
            residue = structure[0][chain_id][res_id]
            is_standard = residue.resname in standard_residues

            # Add original atoms
            for atom in residue:
                element = atom.element if atom.element != " " else atom.name[0]
                original_atoms.append({
                    'element': element,
                    'coord': atom.coord,
                    'residue_name': residue.resname,
                    'residue_num': res_id[1],
                    'atom_name': atom.name,
                    'is_standard': is_standard
                })

            # Add terminal groups for standard residues
            if (chain_id, res_id) in coordinated_standard_residues:
                # Handle N-terminal acetyl group
                n_target, n_bonded = get_local_environment(residue, "N")
                if n_target and n_bonded:
                    result = calculate_terminal_groups_position(n_target, n_bonded)
                    if result and result['type'] == 'acetyl':
                        positions = result['positions']
                        
                        # Add C=O carbon
                        added_groups.append({
                            'element': 'C',
                            'coord': positions['co_carbon'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'CO',
                            'is_standard': is_standard
                        })
                        
                        # Add O
                        added_groups.append({
                            'element': 'O',
                            'coord': positions['oxygen'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'OT',
                            'is_standard': is_standard
                        })
                        
                        # Add methyl carbon
                        added_groups.append({
                            'element': 'C',
                            'coord': positions['methyl_carbon'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'CM',
                            'is_standard': is_standard
                        })
                        
                        # Add methyl hydrogens
                        for i, h_pos in enumerate(positions['hydrogens']):
                            added_groups.append({
                                'element': 'H',
                                'coord': h_pos,
                                'residue_name': residue.resname,
                                'residue_num': res_id[1],
                                'atom_name': f'HM{i+1}',
                                'is_standard': is_standard
                            })
                
                # Handle C-terminal NH2 group
                c_target, c_bonded = get_local_environment(residue, "C")
                if c_target and c_bonded:
                    result = calculate_terminal_groups_position(c_target, c_bonded)
                    if result and result['type'] == 'amino':
                        positions = result['positions']
                        
                        # Add N
                        added_groups.append({
                            'element': 'N',
                            'coord': positions['nitrogen'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'NT',
                            'is_standard': is_standard
                        })
                        
                        # Add H1 and H2
                        added_groups.append({
                            'element': 'H',
                            'coord': positions['hydrogen1'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'HT1',
                            'is_standard': is_standard
                        })
                        
                        added_groups.append({
                            'element': 'H',
                            'coord': positions['hydrogen2'],
                            'residue_name': residue.resname,
                            'residue_num': res_id[1],
                            'atom_name': 'HT2',
                            'is_standard': is_standard
                        })

    # Process all atoms for charge calculations
    for atom in original_atoms:
        atoms_data.append({
            'atomic_number': get_atomic_number(atom['element']),
            'is_standard': 0 if atom['is_standard'] else 0
        })

    for atom in added_groups:
        atoms_data.append({
            'atomic_number': get_atomic_number(atom['element']),
            'is_standard': 1 if atom['is_standard'] else 0
        })

# Write processed charges
    with open('processed_charges.dat', 'w') as charge_file:
        if len(atoms_data) != len(mol2_charges):
            raise Exception(f"Number of atoms mismatch: Structure has {len(atoms_data)}, MOL2 has {len(mol2_charges)}")
            
        for i, atom_data in enumerate(atoms_data):
            charge = 0.000000 if atom_data['is_standard'] == 1 else mol2_charges[i]
            charge_file.write(f"{charge:.6f}\n")

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
        methyl_count = sum(1 for atom in added_groups if atom['atom_name'] in ['CM', 'CN'])
        f.write(f"Total methyl groups added: {methyl_count}\n")
        f.write(f"Total additional atoms added: {len(added_groups)}\n\n")

    return metal_coordination, coordinated_residue_list

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
       'LYS': {'N': -0.3479, 'CA': -0.2400, 'C': 0.7341, 'O': -0.5894},    
       'MET': {'N': -0.4157, 'CA': -0.0237, 'C': 0.5973, 'O': -0.5679},  
       'PHE': {'N': -0.4157, 'CA': -0.0024, 'C': 0.5973, 'O': -0.5679},  
       'PRO': {'N': -0.2548, 'CA': -0.0266, 'C': 0.5896, 'O': -0.5748},  
       'SER': {'N': -0.4157, 'CA': -0.0249, 'C': 0.5973, 'O': -0.5679},  
       'THR': {'N': -0.4157, 'CA': -0.0389, 'C': 0.5973, 'O': -0.5679},  
       'TRP': {'N': -0.4157, 'CA': -0.0275, 'C': 0.5973, 'O': -0.5679},  
       'TYR': {'N': -0.4157, 'CA': -0.0014, 'C': 0.5973, 'O': -0.5679},  
       'VAL': {'N': -0.4157, 'CA': -0.0875, 'C': 0.5973, 'O': -0.5679}
    }
#Read the fixed charges file 
def read_fixed_charges_file(filename):
    fixed_charges_map = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    atom_id = int(parts[0])
                    atom_type = parts[1]
                    residue_name = parts[2]
                    
                    # If the atom type is CA, map it to XC for MOL2 compatibility
                    if atom_type == 'CA':
                        fixed_charges_map[(atom_id, 'XC')] = residue_name
                    else:
                        fixed_charges_map[(atom_id, atom_type)] = residue_name
                        
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Proceeding without fixed charges.")
        return {}
        
    return fixed_charges_map

#Process charges while using reference charges from get_default_reference_charges for N, C, O and CA atoms.
def process_charges(processed_charges_file, mol2_file, target_charge, coordinated_residues, fixed_charges_file):
    # Get reference charges
    reference_charges = get_default_reference_charges()
    
    # Read the fixed charges file
    fixed_charges_map = read_fixed_charges_file(fixed_charges_file)
    
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
                    
                    # Store atom information
                    mol2_atoms.append({
                        'index': atom_index,
                        'id': mol2_atom_id,
                        'type': atom_type,
                        'residue': residue_name
                    })
                    
                    # Check if this atom should have a reference charge
                    key = (mol2_atom_id, atom_type)
                    if key in fixed_charges_map:
                        ref_residue = fixed_charges_map[key]
                        lookup_type = 'CA' if atom_type == 'XC' else atom_type
                        if ref_residue in reference_charges and lookup_type in reference_charges[ref_residue]:
                            new_charges[atom_index] = reference_charges[ref_residue][lookup_type]
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
    parser = argparse.ArgumentParser(description='Analyze metal coordination sites and process charges with methyl group additions')
    parser.add_argument('pdb_file', help='Input PDB file to analyze')
    parser.add_argument('mol2_file', help='Input MOL2 file containing charges')
    parser.add_argument('target_charge', type=float, help='Target total charge for the system')
    parser.add_argument('--metals', nargs='+', 
                       default=['ZN', 'FE', 'CO', 'NI', 'CU', 'MN'],
                       help='List of metals to analyze (default: ZN FE CO NI CU MN)')
    parser.add_argument('--distance', type=float, default=2.5,
                       help='Distance cutoff for coordination in Angstroms (default: 2.5)')

    args = parser.parse_args()

    try:
        metal_sites, coordinated_residues = analyze_and_extract_metal_site(
            args.pdb_file,
            mol2_file=args.mol2_file,
            distance_cutoff=args.distance
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
