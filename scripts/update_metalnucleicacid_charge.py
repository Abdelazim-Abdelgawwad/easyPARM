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

import re
from Bio.PDB import PDBParser, NeighborSearch
import numpy as np
import os

def update_mol2_file(input_mol2, charge_file, output_mol2=None):
    try:
        if output_mol2 is None:
            base, ext = os.path.splitext(input_mol2)
            output_mol2 = f"{os.path.basename(input_mol2)}"

        new_data = {}
        with open(charge_file, 'r') as f:
            for index, line in enumerate(f, 1):
                parts = line.split()
                if len(parts) >= 2:
                    charge = float(parts[0])
                    atom_type = parts[1]
                    new_data[index] = {'charge': charge, 'atom_type': atom_type}

        updated_mol2_lines = []
        with open(input_mol2, 'r') as input_file:
            is_atom_section = False
            for line in input_file:
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    continue
                elif line.startswith("@<TRIPOS>"):
                    is_atom_section = False
                    updated_mol2_lines.append(line)
                    continue

                if is_atom_section and line.strip():
                    atom_match = re.match(
                        r"(\s*\d+\s+)(\S+\s+)(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+\s+)(\d+\s+\S+\s+)(-?\d+\.\d+)",
                        line
                    )
                    if atom_match:
                        atom_id, atom_name, coords, old_atom_type, post_type, old_charge = atom_match.groups()
                        atom_index = int(atom_id.strip())
                        if atom_index in new_data:
                            new_atom_type = new_data[atom_index]['atom_type']
                            new_charge = new_data[atom_index]['charge']
                            updated_line = f"{atom_id}{atom_name}{coords}{new_atom_type:<11}{post_type}{new_charge:8.6f}\n"
                            updated_mol2_lines.append(updated_line)
                        else:
                            updated_mol2_lines.append(line)
                    else:
                        updated_mol2_lines.append(line)
                else:
                    updated_mol2_lines.append(line)

        with open(output_mol2, 'w') as output_file:
            output_file.writelines(updated_mol2_lines)

        return output_mol2

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        raise
    except PermissionError:
        print("Error: Permission denied when accessing files.")
        raise
    except Exception as e:
        print(f"Unexpected error occurred: {e}")
        raise

def batch_update_mol2_files(charges_mapping_file='charges_all.dat'):
    updated_files = []
    try:
        with open(charges_mapping_file, 'r') as mapping_file:
            for line in mapping_file:
                parts = line.strip().split()
                if len(parts) == 2:
                    input_mol2, charge_file = parts
                    try:
                        updated_file = update_mol2_file(input_mol2, charge_file)
                        updated_files.append(updated_file)
                    except Exception as e:
                        print(f"Failed to update {input_mol2}: {e}")

    except FileNotFoundError:
        print(f"Error: Charges mapping file '{charges_mapping_file}' not found.")

    return updated_files

def generate_standard_residue_coordination(input_pdb, output_file="coordinated_residues.txt",
                                           standard_residues=None,
                                           metals=None,
                                           distance_cutoff=2.5,
                                           bond_cutoff=1.75):
    if standard_residues is None:
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

    if metals is None:
        metals = [
            'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Na', 'K', 'Li', 'Rb', 'Cs', 'Mg',
            'Ca', 'Sr', 'Ba', 'V', 'Cr', 'Cd', 'Hg', 'Al', 'Ga', 'In', 'Sn', 'Pb',
            'Bi', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
            'Er', 'Tm', 'Yb', 'Lu', 'Fe2', 'Fe3', 'Fe4', 'Cu1', 'Cu2', 'Mn2', 'Mn3',
            'Mn4', 'Co2', 'Co3', 'Ni2', 'Ni3', 'V2', 'V3', 'V4', 'V5'
        ]

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("nucleic_acid", input_pdb)
    coordination_data = []
    phosphodiester_bond_data = []
    nonstandard_standard_bond_data = []

    metal_coordinated_residues = set()
    standard_residues_to_process = set()

    metal_atoms_found = []
    nonstandard_residues_found = set()
    standard_residues_found = set()

    # Step 1: Find all metal atoms and categorize residues
    metal_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname()
                res_id = residue.get_id()[1]

                if res_name in standard_residues:
                    standard_residues_found.add(f"{res_name}:{res_id}")
                else:
                    nonstandard_residues_found.add(f"{res_name}:{res_id}")

                for atom in residue:
                    atom_name = atom.get_name()

                    is_metal = False
                    for metal in metals:
                        if atom_name == metal or (atom_name.startswith(metal) and atom_name[len(metal):].isdigit()):
                            is_metal = True
                            break

                    if is_metal:
                        metal_atoms.append({
                            'atom': atom,
                            'name': atom_name,
                            'residue': residue,
                            'res_name': res_name,
                            'res_id': res_id,
                            'chain_id': chain.get_id(),
                            'position': atom.coord
                        })
                        metal_atoms_found.append(f"{atom_name} in {res_name}:{res_id}")

    # Step 2: Find standard residues coordinated to metals
    coordination_count = 0
    for metal_info in metal_atoms:
        metal_atom = metal_info['atom']
        metal_name = metal_info['name']
        metal_res_id = metal_info['res_id']
        metal_res_name = metal_info['res_name']
        metal_position = metal_info['position']

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in standard_residues:
                        res_id = residue.get_id()[1]
                        res_name = residue.get_resname()

                        # Don't skip same residue ID if metal is in non-standard residue
                        # Only skip if metal is in a standard residue and it's the same residue
                        if metal_res_name in standard_residues and res_id == metal_res_id:
                            continue

                        for atom in residue:
                            distance = np.linalg.norm(metal_position - atom.coord)

                            if distance <= distance_cutoff:
                                metal_coordinated_residues.add(res_id)
                                standard_residues_to_process.add(res_id)

                                coordination_data.append(
                                    f"bond PRO.{metal_res_id}.{metal_name} "
                                    f"PRO.{res_id}.{atom.get_name()}"
                                )
                                coordination_count += 1

    # Step 3: Find non-standard residues containing metals
    metal_containing_nonstandard = {}

    for metal_info in metal_atoms:
        res_name = metal_info['res_name']
        res_id = metal_info['res_id']

        if res_name not in standard_residues:
            if res_id not in metal_containing_nonstandard:
                metal_containing_nonstandard[res_id] = []
            metal_containing_nonstandard[res_id].append(metal_info)

    # Step 4: Find standard residues that bond with non-standard residues containing metals
    standard_heavy_atoms = []
    standard_atom_info = {}

    nonstandard_heavy_atoms = []
    nonstandard_atom_info = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()

                if res_name in standard_residues:
                    for atom in residue:
                        element = atom.element if hasattr(atom, 'element') else atom.get_name()[0]
                        if element != 'H':
                            standard_heavy_atoms.append(atom)
                            standard_atom_info[atom] = (res_id, atom.get_name())

                elif res_id in metal_containing_nonstandard:
                    for atom in residue:
                        is_metal_atom = False
                        for metal_info in metal_containing_nonstandard[res_id]:
                            if atom == metal_info['atom']:
                                is_metal_atom = True
                                break

                        if not is_metal_atom:
                            element = atom.element if hasattr(atom, 'element') else atom.get_name()[0]
                            if element != 'H':
                                nonstandard_heavy_atoms.append(atom)
                                nonstandard_atom_info[atom] = (res_id, atom.get_name())

    # Find bonds between non-standard residues with metals and standard residues
    if nonstandard_heavy_atoms and standard_heavy_atoms:
        ns = NeighborSearch(standard_heavy_atoms)

        nonstandard_residues_linked_to_standard = set()

        for atom in nonstandard_heavy_atoms:
            nonstandard_res_id, atom_name = nonstandard_atom_info[atom]

            nearby_atoms = ns.search(atom.coord, bond_cutoff, level='A')

            for nearby_atom in nearby_atoms:
                if nearby_atom in standard_atom_info:
                    standard_res_id, standard_atom_name = standard_atom_info[nearby_atom]

                    nonstandard_standard_bond_data.append(
                        f"bond PRO.{nonstandard_res_id}.{atom_name} "
                        f"PRO.{standard_res_id}.{standard_atom_name}"
                    )

                    standard_residues_to_process.add(standard_res_id)
                    nonstandard_residues_linked_to_standard.add(nonstandard_res_id)
    # Step 5: Find phosphodiester bonds in nucleic acid backbone
    backbone_atoms = []
    backbone_atom_info = {}

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()

                if res_name in standard_residues:
                    for atom in residue:
                        atom_name = atom.get_name()
                        if atom_name in ["P", "O3'"]:
                            backbone_atoms.append(atom)
                            backbone_atom_info[atom] = (res_id, atom_name, chain_id)

    phosphodiester_bonds_set = set()

    if backbone_atoms:
        ns_backbone = NeighborSearch(backbone_atoms)

        for atom in backbone_atoms:
            res_id, atom_name, chain_id = backbone_atom_info[atom]

            if res_id in standard_residues_to_process:
                nearby_atoms = ns_backbone.search(atom.coord, bond_cutoff, level='A')

                for nearby_atom in nearby_atoms:
                    if nearby_atom in backbone_atom_info:
                        nearby_res_id, nearby_atom_name, nearby_chain_id = backbone_atom_info[nearby_atom]

                        if nearby_res_id != res_id and chain_id == nearby_chain_id:
                            if (atom_name == "O3'" and nearby_atom_name == "P") or \
                               (atom_name == "P" and nearby_atom_name == "O3'"):
                                
                                if atom_name == "O3'" and nearby_atom_name == "P":
                                    if res_id < nearby_res_id:
                                        bond = f"bond PRO.{res_id}.O3' PRO.{nearby_res_id}.P"
                                    else:
                                        bond = f"bond PRO.{nearby_res_id}.P PRO.{res_id}.O3'"
                                else:
                                    if nearby_res_id < res_id:
                                        bond = f"bond PRO.{nearby_res_id}.O3' PRO.{res_id}.P"
                                    else:
                                        bond = f"bond PRO.{res_id}.P PRO.{nearby_res_id}.O3'"

                                if res_id in standard_residues_to_process or nearby_res_id in standard_residues_to_process:
                                    phosphodiester_bonds_set.add(bond)

    phosphodiester_bond_data = [f"{bond}" for bond in sorted(phosphodiester_bonds_set)]
    
    # Combine all data with proper section headers
    all_data = []

    if coordination_data:
        if not all_data:  # Add header only once
            all_data.append("# All bond information for tLeap input file")
        all_data.append("# Metal-nucleotide coordination bonds:")
        all_data.extend(coordination_data)
        all_data.append("")

    if nonstandard_standard_bond_data:
        if not all_data:  # Add header only once
            all_data.append("# All bond information for tLeap input file")
        all_data.append("# Ligand-nucleotide covalent bonds:")
        all_data.extend(nonstandard_standard_bond_data)
        all_data.append("")

    if phosphodiester_bond_data:
        if not all_data:  # Add header only once
            all_data.append("# All bond information for tLeap input file")
        all_data.append("# Phosphodiester backbone bonds:")
        all_data.extend(phosphodiester_bond_data)

    # Write output if any data was found
    if all_data:
        with open(output_file, "w") as f:
            f.write("\n".join(all_data))
    else:
        print("\n=== No coordination data found ===")

if __name__ == "__main__":
    input_pdb = "easyPARM.pdb"

    batch_update_mol2_files('charges_all.dat')
    generate_standard_residue_coordination(input_pdb, output_file="coordinated_residues.txt", distance_cutoff=2.5)
