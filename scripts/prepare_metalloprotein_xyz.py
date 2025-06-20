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
#                              |  $$$$$$/              Ver. 4.00 - 8 June 2025                                    #
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

# Check the surronding atom to better adding the capping atom
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

# To reverse the adding part
def rotate_vector(v, axis, angle):
    cos_ang = np.cos(angle)
    sin_ang = np.sin(angle)
    return (v * cos_ang + 
            np.cross(axis, v) * sin_ang + 
            axis * np.dot(axis, v) * (1 - cos_ang))

# Add hydrogen to methyl group
def calculate_methyl_hydrogens(c_pos, target_coord, methyl_direction, h_bond_length=1.09):
    h_bond_angle = np.radians(109.5)
    perp1 = get_perpendicular_vector(methyl_direction)
    perp2 = np.cross(methyl_direction, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)
    
    h_positions = []
    for i in range(3):
        rotation = i * np.radians(120)
        h_direction = (np.cos(h_bond_angle) * -methyl_direction + 
                      np.sin(h_bond_angle) * (np.cos(rotation) * perp1 + 
                                            np.sin(rotation) * perp2))
        h_pos = c_pos + h_direction * h_bond_length
        h_positions.append(h_pos)
    return h_positions

#Calculate vectors defining the amide plane based on existing bonds
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

#Calculate optimal acetyl group position based on local bonding environment
#with improved amide geometry and steric considerations
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

#Calculate acetyl group positions with improved geometric constraints
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


def analyze_and_extract_metal_site(input_pdb, metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5'], distance_cutoff=2.6):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    standard_residues = {
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
        "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
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
        bond_cutoff = 1.9  # Covalent bond distance cutoff in Å
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
    # Create a subset of atoms only from residues we're extracting
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

        # Calculate terminal modifications if needed
        if res_key in coordinated_residues:
            # Check if N atom needs acetylation
            n_atom, n_bonded = get_local_environment(residue, "N", max_bond_distance=1.6)
            if n_atom is not None:
                needs_acetyl = True
                
                # Check if N already has cross-residue bonds within the extracted residues
                if "N" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["N"]:
                    # N has bonds to another residue within our selection
                    needs_acetyl = False
                
                # Alternatively, check if N already has 3+ total bonds (already capped or part of chain)
                if len(n_bonded) >= 3:
                    needs_acetyl = False
                    
                if needs_acetyl:
                    c_pos, o_pos, ch3_pos, h_positions = calculate_acetyl_group(n_atom, n_bonded)
                    if c_pos is not None:
                        terminal_modifications.extend([
                            {
                                'element': 'C', 'coord': c_pos, 'name': 'CAC',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'O', 'coord': o_pos, 'name': 'OAC',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            },
                            {
                                'element': 'C', 'coord': ch3_pos, 'name': 'CME',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            }
                        ])
                        for i, h_pos in enumerate(h_positions):
                            terminal_modifications.append({
                                'element': 'H', 'coord': h_pos, 'name': f'HM{i+1}',
                                'resname': residue.resname, 'resnum': residue.get_id()[1],
                                'is_capping': True
                            })

            # Check if C atom needs amination
            c_atom, c_bonded = get_local_environment(residue, "C", max_bond_distance=1.6)
            if c_atom is not None:
                needs_amination = True
                
                # Check if C already has cross-residue bonds within the extracted residues
                if "C" in cross_residue_bonds[res_key] and cross_residue_bonds[res_key]["C"]:
                    # C has bonds to another residue within our selection
                    needs_amination = False
                
                # For C atom in standard residues, it should have 3 bonds if capped (CA, O, NH2)
                # or 3 bonds if part of a peptide chain (CA, O, next residue's N)
                # If it only has 2 bonds (likely CA and O), it needs capping
                if len(c_bonded) >= 3:
                    needs_amination = False
                    
                if needs_amination:
                    n_pos, h_positions = calculate_amino_group(c_atom, c_bonded)
                    if n_pos is not None:
                        terminal_modifications.append({
                            'element': 'N', 'coord': n_pos, 'name': 'NT',
                            'resname': residue.resname, 'resnum': residue.get_id()[1],
                            'is_capping': True
                        })
                        for i, h_pos in enumerate(h_positions):
                            terminal_modifications.append({
                                'element': 'H', 'coord': h_pos, 'name': f'HT{i+1}',
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
    parser = argparse.ArgumentParser(description='Analyze metal coordination sites and add terminal modifications')
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
        sys.exit(1)

if __name__ == "__main__":
    main()

