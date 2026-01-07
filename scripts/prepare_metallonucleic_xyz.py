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
#                              |  $$$$$$/              Ver. 4.15 - 17 October 2025                                #
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

# Check the surrounding atom to better add the capping atom
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

#Calculate capping group for O3' terminal with proper tetrahedral geometry:
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
    
    # Create proper tetrahedral geometry around phosphorus
    # The key is to use the tetrahedral angle (109.47°) from the P-O3' bond axis
    tet_angle = np.radians(109.47)
    
    # Get two perpendicular vectors to avg_direction
    perp1 = get_perpendicular_vector(avg_direction)
    perp2 = np.cross(avg_direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    # Double-bonded oxygen (P=O): place opposite to O3' at tetrahedral angle
    o_double_bond = 1.48  # P=O bond length
    # Rotate away from O3' by tetrahedral angle
    o_double_direction = rotate_vector(-avg_direction, perp1, 0)  # Directly opposite
    o_double_pos = p_pos + o_double_direction * o_double_bond
    
    # Two single-bonded oxygens: place at 120° apart around the P-O3' axis
    # and at tetrahedral angle from the P-O3' bond
    o_single_bond = 1.61  # P-O single bond
    
    # Create a cone of directions at tetrahedral angle from -avg_direction
    # First single-bonded O: rotate by tetrahedral angle, then 0° around axis
    base_direction = rotate_vector(-avg_direction, perp1, tet_angle)
    o1_pos = p_pos + base_direction * o_single_bond
    
    # Second single-bonded O: rotate around the P-O3' axis by 120°
    o2_direction = rotate_vector(base_direction, -avg_direction, np.radians(120))
    o2_pos = p_pos + o2_direction * o_single_bond
    
    # Third position would be at 240° for the double-bonded oxygen
    # Let's recalculate to ensure proper spacing
    o_double_direction_corrected = rotate_vector(base_direction, -avg_direction, np.radians(240))
    o_double_pos = p_pos + o_double_direction_corrected * o_double_bond
    
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


def analyze_and_extract_metal_site(input_pdb, metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5'], distance_cutoff=2.6):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('nucleic_acid', input_pdb)

    # Standard nucleic acid residues (DNA and RNA)
    standard_residues = {
        # DNA
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

    metal_sites = {}
    residues_to_extract = set()
    coordinated_residues = set()
    all_atoms = []
    terminal_modifications = []

    # First pass: Find metal sites and coordinating residues
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element in metals:
                        metal_coord = atom.coord
                        metal_key = f"{atom.element}_{chain.id}_{residue.get_id()[1]}"
                        metal_sites[metal_key] = {'metal': atom, 'coordinating': []}

                        for chain2 in model:
                            for residue2 in chain2:
                                is_coordinating = False
                                for atom2 in residue2:
                                    distance = np.linalg.norm(metal_coord - atom2.coord)
                                    if distance <= distance_cutoff:
                                        is_coordinating = True
                                        metal_sites[metal_key]['coordinating'].append(residue2)
                                        residues_to_extract.add((chain2.id, residue2.get_id()))

                                if is_coordinating and residue2.get_resname() in standard_residues:
                                    coordinated_residues.add((chain2.id, residue2.get_id()))

        # Find standard residues linked to non-standard coordinating residues
        heavy_atoms = [atom for chain in model for residue in chain for atom in residue if atom.element != 'H']
        ns = PDB.NeighborSearch(heavy_atoms)
        bond_cutoff = 1.9  # Covalent bond distance cutoff in Ã…
        
        for metal_key in metal_sites:
            for coord_res in metal_sites[metal_key]['coordinating']:
                if coord_res.resname not in standard_residues:  # Non-standard residue
                    coord_heavy_atoms = [atom for atom in coord_res if atom.element != 'H']
                    linked_residues = set()
                    for atom in coord_heavy_atoms:
                        nearby_atoms = ns.search(atom.coord, bond_cutoff, level='A')
                        for nearby_atom in nearby_atoms:
                            nearby_res = nearby_atom.get_parent()
                            if (nearby_res != coord_res and
                                nearby_res.resname in standard_residues):
                                linked_residues.add(nearby_res)
                    for linked_res in linked_residues:
                        res_key = (linked_res.get_parent().id, linked_res.get_id())
                        if res_key not in residues_to_extract:
                            residues_to_extract.add(res_key)
                            coordinated_residues.add(res_key)

    # Improved cross-residue bond detection focused on the extracted residues
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
    bond_cutoff = 1.9  # Covalent bond distance cutoff
    
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

    # Second pass: Extract atoms and calculate modifications
    for res_key, residue in extracted_residue_objects.items():
        # Add original atoms
        for atom in residue:
            all_atoms.append({
                'element': atom.element if atom.element != " " else atom.name[0],
                'coord': atom.coord,
                'name': atom.name,
                'resname': residue.resname,
                'resnum': residue.get_id()[1],
                'is_capping': False
            })

        # Calculate terminal modifications if needed (only for standard nucleic acid residues)
        if res_key in coordinated_residues:
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
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'O', 'coord': cap_data['O_double'], 'name': 'O1C',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'O', 'coord': cap_data['O_single1'], 'name': 'O2C',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'O', 'coord': cap_data['O_single2'], 'name': 'O3C',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'C', 'coord': cap_data['CH3'], 'name': 'CMC',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            }
                        ])
                        for i, h_pos in enumerate(cap_data['H_positions']):
                            terminal_modifications.append({
                                'element': 'H', 'coord': h_pos, 'name': f'HMC{i+1}',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
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
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'C', 'coord': cap_data['CH3'], 'name': 'CPC',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            }
                        ])
                        for i, h_pos in enumerate(cap_data['H_positions']):
                            terminal_modifications.append({
                                'element': 'H', 'coord': h_pos, 'name': f'HPC{i+1}',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            })

    # Combine original atoms and modifications
    all_atoms.extend(terminal_modifications)

    # Write XYZ file with conditional formatting
    with open('initial_structure.xyz', 'w') as f:
        f.write(f"{len(all_atoms)}\n\n")
        for atom in all_atoms:
            if atom['is_capping']:
                f.write(f"{atom['element']:2} {atom['coord'][0]:10.6f} {atom['coord'][1]:10.6f} {atom['coord'][2]:10.6f} {atom['name']:4} {atom['resname']:3} {atom['resnum']:4}  CAPPING\n")
            else:
                f.write(f"{atom['element']:2} {atom['coord'][0]:10.6f} {atom['coord'][1]:10.6f} {atom['coord'][2]:10.6f} {atom['name']:4} {atom['resname']:3} {atom['resnum']:4}\n")

    return metal_sites

def main():
    parser = argparse.ArgumentParser(description='Analyze metal coordination sites in nucleic acids and add terminal modifications')
    parser.add_argument('pdb_file', help='Input PDB file to analyze')

    args = parser.parse_args()

    try:
        metal_sites = analyze_and_extract_metal_site(
            args.pdb_file
        )
        
    except FileNotFoundError:
        print(f"Error: File '{args.pdb_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
