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
#                              |  $$$$$$/              Ver. 4.20 - 1 January 2026                                 #
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

# Calculate methyl group position for phosphate capping
def calculate_methyl_hydrogens(c_pos, direction_vector, h_bond_length=1.09):
    h_bond_angle = np.radians(109.5)
    perp1 = get_perpendicular_vector(direction_vector)
    perp2 = np.cross(direction_vector, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    h_positions = []
    for i in range(3):
        rotation = i * np.radians(120)
        h_direction = (np.cos(h_bond_angle) * direction_vector + 
                      np.sin(h_bond_angle) * (np.cos(rotation) * perp1 + 
                                            np.sin(rotation) * perp2))
        h_pos = c_pos + h_direction * h_bond_length
        h_positions.append(h_pos)
    return h_positions

#Calculate capping group for O3' terminal:
def calculate_phosphate_cap_for_o3prime(o3_atom, bonded_atoms):
    """
              O
              ||
    O3'---P---O---CH3
              |
              O
    """
    if not bonded_atoms:
        return None
        
    o3_coord = o3_atom.coord
    
    # Calculate average direction from bonded atoms
    avg_direction = np.zeros(3)
    for atom in bonded_atoms:
        direction = atom.coord - o3_coord
        direction = direction / np.linalg.norm(direction)
        avg_direction += direction
    
    avg_direction = -avg_direction / len(bonded_atoms)
    avg_direction = avg_direction / np.linalg.norm(avg_direction)
    
    # Place phosphorus
    p_o_bond = 1.61  # P-O bond length
    p_pos = o3_coord + avg_direction * p_o_bond
    
    # Calculate three oxygen positions around phosphorus (tetrahedral)
    # One double-bonded O, two single-bonded O (one will connect to CH3)
    perp1 = get_perpendicular_vector(avg_direction)
    perp2 = np.cross(avg_direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # Double-bonded oxygen (opposite to O3')
    o_double_bond = 1.48  # P=O bond length
    o_double_pos = p_pos + (-avg_direction) * o_double_bond
    
    # Two single-bonded oxygens at tetrahedral angles
    tet_angle = np.radians(109.5)
    o_single_bond = 1.61  # P-O single bond
    
    # First single-bonded O
    o1_direction = rotate_vector(-avg_direction, perp1, tet_angle)
    o1_pos = p_pos + o1_direction * o_single_bond
    
    # Second single-bonded O (this will connect to CH3)
    o2_direction = rotate_vector(-avg_direction, perp1, -tet_angle)
    o2_pos = p_pos + o2_direction * o_single_bond
    
    # Calculate methyl group attached to o2
    c_o_bond = 1.43  # C-O bond length
    c_direction = o2_direction / np.linalg.norm(o2_direction)
    ch3_pos = o2_pos + c_direction * c_o_bond
    
    # Calculate hydrogen positions for methyl group
    h_positions = calculate_methyl_hydrogens(ch3_pos, c_direction)
    
    return {
        'P': p_pos,
        'O_double': o_double_pos,
        'O_single1': o1_pos,
        'O_single2': o2_pos,
        'CH3': ch3_pos,
        'H_positions': h_positions
    }

#Calculate capping group for 5' phosphate terminal:
def calculate_methyl_cap_for_phosphate(p_atom, bonded_atoms):
    """
    P---O---CH3
    """
    if not bonded_atoms:
        return None
        
    p_coord = p_atom.coord
    
    # Calculate direction away from bonded atoms
    avg_direction = np.zeros(3)
    for atom in bonded_atoms:
        direction = atom.coord - p_coord
        direction = direction / np.linalg.norm(direction)
        avg_direction += direction
    
    avg_direction = -avg_direction / len(bonded_atoms)
    avg_direction = avg_direction / np.linalg.norm(avg_direction)
    
    # Place oxygen
    p_o_bond = 1.61  # P-O bond length
    o_pos = p_coord + avg_direction * p_o_bond
    
    # Place methyl carbon
    c_o_bond = 1.43  # C-O bond length
    c_pos = o_pos + avg_direction * c_o_bond
    
    # Calculate hydrogen positions
    h_positions = calculate_methyl_hydrogens(c_pos, avg_direction)
    
    return {
        'O': o_pos,
        'CH3': c_pos,
        'H_positions': h_positions
    }

def analyze_and_extract_metal_site(input_pdb, mol2_file, metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5'], distance_cutoff=2.6):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('nucleic_acid', input_pdb)
    
    if isinstance(metals, list):
        metals = set(metals)
    elif not isinstance(metals, set):
        metals = {'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 
                'NA', 'K', 'CA', 'LI', 'RB', 'CS', 'MG', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 
                'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}
    
    # Standard nucleic acid residues (DNA and RNA)
    standard_residues = {
        "DA", "DT", "DC", "DG", "DU",
        # RNA
        "A", "U", "C", "G",
        # Additional AMBER nucleic acid residues
        "RA", "RU", "RC", "RG",  # RNA
        "DAN", "DTN", "DCN", "DGN",  # Deoxy forms
        "A3", "A5", "AN", "C3", "C5", "CN",  # 3' and 5' terminal
        "G3", "G5", "GN", "U3", "U5", "UN",
        "DA3", "DA5", "DAN", "DC3", "DC5", "DCN",
        "DG3", "DG5", "DGN", "DT3", "DT5", "DTN",
        "OHE", "ADE", "GUA", "CYT", "THY", "URA",
        "RA3", "RA5", "RU3", "RU5", "RG3", "RG5", "RC3", "RC5"
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
    bond_cutoff = 1.9  # Covalent bond distance cutoff in Ã…
    
    # Find standard residues linked to non-standard coordinating residues
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
                        if atom.element != 'H':  # Only include heavy atoms
                            extracted_atoms.append(atom)
    
    # Use NeighborSearch only on atoms within the extracted residues
    ns_extracted = PDB.NeighborSearch(extracted_atoms)
    
    # Initialize bond dictionary for extracted residues
    cross_residue_bonds = {}
    for res_key in residues_to_extract:
        cross_residue_bonds[res_key] = {}
    
    # Detect all cross-residue bonds within the extracted residues
    for atom in extracted_atoms:
        res_key = (atom.get_parent().get_parent().id, atom.get_parent().get_id())
        atom_name = atom.name
        
        # Initialize entry for this atom if not exists
        if atom_name not in cross_residue_bonds[res_key]:
            cross_residue_bonds[res_key][atom_name] = []
        
        # Look for nearby atoms that could form bonds
        nearby_atoms = ns_extracted.search(atom.coord, bond_cutoff, level='A')
        for nearby_atom in nearby_atoms:
            if nearby_atom == atom:
                continue  # Skip self
                
            nearby_res_key = (nearby_atom.get_parent().get_parent().id, 
                             nearby_atom.get_parent().get_id())
            
            # Only record if it's a different residue
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
                    'is_standard': 0 if is_standard else 0
                })
            
            # Add terminal groups for standard nucleic acid residues that coordinate with metals
            if (chain_id, res_id) in coordinated_standard_residues:
                res_key = (chain_id, res_id)
                
                # Check if O3' needs phosphate capping
                o3_atom, o3_bonded = get_local_environment(residue, "O3'", max_bond_distance=1.7)
                if o3_atom is not None:
                    needs_o3_cap = True
                    
                    # Check if O3' already has cross-residue bonds (connected to next residue's P)
                    if "O3'" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["O3'"]:
                        needs_o3_cap = False
                    
                    # Check if O3' already has 2+ bonds (already bonded to P)
                    if len(o3_bonded) >= 2:
                        needs_o3_cap = False
                        
                    if needs_o3_cap:
                        cap_data = calculate_phosphate_cap_for_o3prime(o3_atom, o3_bonded)
                        if cap_data is not None:
                            terminal_modifications.extend([
                                {
                                    'element': 'P', 'coord': cap_data['P'], 'name': 'PC',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                },
                                {
                                    'element': 'O', 'coord': cap_data['O_double'], 'name': 'O1C',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                },
                                {
                                    'element': 'O', 'coord': cap_data['O_single1'], 'name': 'O2C',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                },
                                {
                                    'element': 'O', 'coord': cap_data['O_single2'], 'name': 'O3C',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                },
                                {
                                    'element': 'C', 'coord': cap_data['CH3'], 'name': 'CMC',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                }
                            ])
                            for i, h_pos in enumerate(cap_data['H_positions']):
                                terminal_modifications.append({
                                    'element': 'H', 'coord': h_pos, 'name': f'HMC{i+1}',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                })

                # Check if 5' phosphate (P) needs methyl capping
                p_atom, p_bonded = get_local_environment(residue, "P", max_bond_distance=1.7)
                if p_atom is not None:
                    needs_p_cap = True
                    
                    # Check if P already has cross-residue bonds (connected to previous residue's O3')
                    if "P" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["P"]:
                        # Check if any of the bonded atoms is an O3' from another residue
                        for bonded_atom in cross_residue_bonds[res_key]["P"]:
                            if bonded_atom.name == "O3'":
                                needs_p_cap = False
                                break
                    
                    # P typically has 4 bonds: 3 oxygens + 1 to O5'
                    # If it has all expected bonds within the residue + cross-residue, no cap needed
                    # Count oxygen bonds
                    o_bonds = sum(1 for atom in p_bonded if atom.element == 'O')
                    if o_bonds >= 4:  # Full coordination
                        needs_p_cap = False
                        
                    if needs_p_cap:
                        cap_data = calculate_methyl_cap_for_phosphate(p_atom, p_bonded)
                        if cap_data is not None:
                            terminal_modifications.extend([
                                {
                                    'element': 'O', 'coord': cap_data['O'], 'name': 'OPC',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                },
                                {
                                    'element': 'C', 'coord': cap_data['CH3'], 'name': 'CPC',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                }
                            ])
                            for i, h_pos in enumerate(cap_data['H_positions']):
                                terminal_modifications.append({
                                    'element': 'H', 'coord': h_pos, 'name': f'HPC{i+1}',
                                    'resname': residue.resname, 'resnum': res_id[1],
                                    'is_capping': True
                                })
    
    # Add terminal modifications to all_atoms
    all_atoms.extend(terminal_modifications)
    
    # Add terminal modifications to atoms_data for charge calculations
    for mod in terminal_modifications:
        atoms_data.append({
            'atomic_number': get_atomic_number(mod['element']),
            'is_standard': 1  # Capping groups always marked as standard
        })
    
    # Calculate the total charge contribution from capping atoms
    # For phosphate cap: -1 charge (P=O^- and two O^-)
    # For methyl cap: 0 charge (neutral methyl group)
    total_cap_charge = 0.0
    phosphate_cap_count = sum(1 for atom in terminal_modifications if atom['name'] == 'PC')  # Count phosphate caps
    total_cap_charge = -phosphate_cap_count  # Each phosphate cap contributes -1 charge

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
                       f"at {res['distance']} Ã…\n")
            f.write("\n")
        
        f.write("\nCAPPING Group Additions:\n")
        f.write("=====================\n")
        phosphate_cap_count = sum(1 for atom in terminal_modifications if atom['name'] == 'PC')
        methyl_cap_count = sum(1 for atom in terminal_modifications if atom['name'] == 'OPC')
        f.write(f"Total O3' phosphate caps added: {phosphate_cap_count}\n")
        f.write(f"Total 5' phosphate methyl caps added: {methyl_cap_count}\n")
        f.write(f"Total additional atoms added: {len(terminal_modifications)}\n\n")
    
    return metal_coordination, coordinated_residue_list, total_cap_charge


#Provides default charge for N,C,O and CA atoms from amber force field.
def get_default_reference_charges():
    return {
    'A': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'A3': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'A5': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'AN': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'C': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'C3': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'C5': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'CN': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'G': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'G3': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'G5': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'GN': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'OHE': {'HOP3': 0.3129, 'OP3': -0.621},
    'U': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'RA': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RA3': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.6541},
    'RA5': {"O5'": -0.6223, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RC': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RC3': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.6541},
    'RC5': {"O5'": -0.6223, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RG': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RG3': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.6541},
    'RG5': {"O5'": -0.6223, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RU': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'RU3': {'P': 1.1662, "O5'": -0.4989, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.6541},
    'RU5': {"O5'": -0.6223, "C5'": 0.0558, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "O3'": -0.5246},
    'U3': {'P': 1.1662, 'OP1': -0.776, 'OP2': -0.776, "O5'": -0.4989, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'U5': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.5246},
    'UN': {"O5'": -0.6223, "C5'": 0.0558, "H5'": 0.0679, "H5''": 0.0679, "C4'": 0.1065, "H4'": 0.1174, "O4'": -0.3548, "C3'": 0.2022, "H3'": 0.0615, "C2'": 0.067, "H2'": 0.0972, "O3'": -0.6541},
    'DA': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DA3': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DA5': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DAN': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DC': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DC3': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DC5': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DCN': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DG': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DG3': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DG5': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DGN': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DT': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DT3': {'P': 1.1659, 'OP1': -0.7761, 'OP2': -0.7761, "O5'": -0.4954, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549},
    'DT5': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.5232},
    'DTN': {"O5'": -0.6318, "C5'": -0.0069, "H5'": 0.0754, "H5''": 0.0754, "C4'": 0.1629, "H4'": 0.1176, "O4'": -0.3691, "C3'": 0.0713, "H3'": 0.0985, "C2'": -0.0854, "H2'": 0.0718, "H2''": 0.0718, "O3'": -0.6549}}

#Read the fixed charges file 
def read_fixed_charges_file(filename):
    fixed_charges_map = {}

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()

                # Expect: atom_id, atom_name_ref, atom_name_out, residue_name
                if len(parts) >= 4:
                    atom_id = int(parts[0])
                    atom_name_ref = parts[1]    # Reference name for lookup (e.g., "O3'")
                    atom_name_out = parts[2]    # Name in MOL2 file (e.g., "O2")
                    residue_name = parts[3]

                    # Map (atom_id, atom_name_out) -> (atom_name_ref, residue_name)
                    # This allows us to look up using the MOL2 name and get the reference name
                    fixed_charges_map[(atom_id, atom_name_out)] = (atom_name_ref, residue_name)

    except FileNotFoundError:
        print(f"Warning: {filename} not found. Proceeding without fixed charges.")
        return {}

    return fixed_charges_map

#Process charges while using reference charges from get_default_reference_charges for N, C, O and CA atoms.
def process_charges(processed_charges_file, mol2_file, target_charge, coordinated_residues, fixed_charges_file, total_cap_charge=0.0):

    # Get reference charges
    reference_charges = get_default_reference_charges()
    
    # Read the fixed charges file
    # Now returns: (atom_id, atom_name_out) -> (atom_name_ref, residue_name)
    fixed_charges_map = read_fixed_charges_file(fixed_charges_file)
    target_charge = target_charge - total_cap_charge 
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
                    atom_type = parts[5]       # This is atom_name_out in MOL2
                    residue_name = parts[7]
                    
                    # Store atom information
                    mol2_atoms.append({
                        'index': atom_index,
                        'id': mol2_atom_id,
                        'type': atom_type,
                        'residue': residue_name
                    })
                    
                    # Check if this atom should have a reference charge
                    # Look up using (mol2_atom_id, atom_type) where atom_type is from MOL2
                    key = (mol2_atom_id, atom_type)
                    
                    if key in fixed_charges_map:
                        # Get the reference name and residue for charge lookup
                        atom_name_ref, ref_residue = fixed_charges_map[key]
                        
                        # Now look up the charge using atom_name_ref
                        if ref_residue in reference_charges and atom_name_ref in reference_charges[ref_residue]:
                            new_charges[atom_index] = reference_charges[ref_residue][atom_name_ref]
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
        metal_sites, coordinated_residues, total_cap_charge = analyze_and_extract_metal_site(
            args.pdb_file,
            mol2_file=args.mol2_file,
        )
        
        process_charges('processed_charges.dat', args.mol2_file, args.target_charge, coordinated_residues, fixed_charges_file="fixed_charges.dat", total_cap_charge=total_cap_charge) 
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
