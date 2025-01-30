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
import os 

def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        atom_types = {}
        is_atom_section = False
        
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_type = parts[5]
                    atom_types[atom_id] = atom_type
        
        return atom_types
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

#Read atom type mappings from library file
def read_aminofb15_lib(file_path):
    residue_atom_types = {}
    current_residue = None
    
    try:
        with open(file_path, 'r') as file:
            reading_atoms = False
            for line in file:
                if '!entry.' in line and '.unit.atoms table' in line:
                    current_residue = line.split('.')[1]
                    reading_atoms = True
                    residue_atom_types[current_residue] = {}
                elif '!entry.' in line and '.unit.atomspertinfo table' in line:
                    reading_atoms = False
                elif reading_atoms and line.strip() and '"' in line:
                    parts = line.strip().split('"')
                    if len(parts) >= 4:
                        orig_atom_name = parts[1]
                        new_atom_type = parts[3]
                        residue_atom_types[current_residue][orig_atom_name] = new_atom_type
        
        return residue_atom_types
    except Exception as e:
        print(f"Error reading library file: {e}")
        return {}

#Write initial PDB with MOL2 atom types
def write_pdb_with_mol2_types(structure, output_file, atom_types):
    atom_counter = 1
    with open(output_file, 'w') as f:
        for model in structure:
            for chain in model:
                for residue in chain:
                    hetero, resnum, icode = residue.get_id()
                    for atom in residue:
                        record_type = "HETATM" if residue.get_id()[0] != " " else "ATOM"
                        mol2_type = atom_types.get(atom_counter, atom.name)
                        x, y, z = atom.get_coord()
                        resname = residue.get_resname()
                        chainid = chain.get_id()
                        occupancy = atom.get_occupancy()
                        bfactor = atom.get_bfactor()
                        element = atom.element if hasattr(atom, 'element') else '  '
                        
                        line = f"{record_type:<6}{atom_counter:>5}  {mol2_type:<3} {resname:>3} {chainid}{resnum:>4}    "
                        line += f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}          {element:>2}"
                        f.write(line + '\n')
                        atom_counter += 1

#Update atom types in QM.pdb based on REFQM.pdb and library
def update_atom_types_from_library(qm_pdb, ref_pdb, lib_types):
    with open(ref_pdb, 'r') as ref_file:
        ref_lines = ref_file.readlines()
    
    updated_lines = []
    with open(qm_pdb, 'r') as qm_file:
        qm_lines = qm_file.readlines()
        
        for qm_line, ref_line in zip(qm_lines, ref_lines):
            if qm_line.startswith(('ATOM', 'HETATM')):
                # Get residue name and original atom name from reference PDB
                ref_resname = ref_line[17:20].strip()
                ref_atomname = ref_line[12:16].strip()
                
                # Look up new atom type in library
                if ref_resname in lib_types and ref_atomname in lib_types[ref_resname]:
                    new_type = lib_types[ref_resname][ref_atomname]
                    # Replace atom type in QM.pdb line (columns 13-16)
                    new_line = qm_line[:12] + f"{new_type:<3}" + qm_line[16:]
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(qm_line)
            else:
                updated_lines.append(qm_line)
    
    # Write updated QM.pdb
    with open(qm_pdb, 'w') as f:
        f.writelines(updated_lines)

#Update atom types in MOL2 file using atom types from QM.pdb while preserving the original formatting.
def update_mol2_with_qm_types(mol2_file, qm_pdb, output_mol2):

    # Read QM.pdb atom types
    qm_atom_types = {}
    with open(qm_pdb, 'r') as qm_file:
        for line in qm_file:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    atom_number = int(line[6:11].strip())
                    atom_type = line[12:16].strip()
                    qm_atom_types[atom_number] = atom_type
                except ValueError:
                    print(f"Warning: Skipping malformed QM.pdb line: {line.strip()}")

    # Read and update MOL2 file
    updated_mol2_lines = []
    try:
        with open(mol2_file, 'r') as mol2:
            lines = mol2.readlines()
            
            is_atom_section = False
            for line in lines:
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    continue
                
                if line.startswith("@<TRIPOS>BOND"):
                    is_atom_section = False
                
                if is_atom_section and line.strip():
                    try:
                        # Use regex to parse MOL2 atom line while preserving spacing
                        atom_match = re.match(
                            r"(\s*\d+\s+)(\S+)(\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+)(\s+\d+\s+\S+\s+-?\d+\.\d+)", 
                            line
                        )
                        
                        if atom_match:
                            # Extract groups
                            pre_id, atom_name, coords, atom_type, post_type = atom_match.groups()
                            atom_id = int(pre_id.strip())

                            # Update atom type if found in QM.pdb
                            if atom_id in qm_atom_types:
                                atom_type = f"{qm_atom_types[atom_id]:<2}"  # Left-align to match MOL2 format

                            updated_line = f"{pre_id}{atom_name}{coords}{atom_type}{post_type}\n"
                            updated_mol2_lines.append(updated_line)
                        else:
                            updated_mol2_lines.append(line)
                    except Exception as e:
                        print(f"Warning: Skipping malformed MOL2 line: {line.strip()} ({e})")
                        updated_mol2_lines.append(line)  # Keep the original line
                else:
                    updated_mol2_lines.append(line)
        
        # Write updated MOL2 file
        with open(output_mol2, 'w') as output_file:
            output_file.writelines(updated_mol2_lines)
         
    except Exception as e:
        print(f"Error updating MOL2 file: {e}")

def analyze_and_extract_metal_site(input_pdb, mol2_file="NEW_COMPLEX.mol2", lib_file="amber_refernce.lib", 
                                   output_pdb="QM.pdb", distance_cutoff=2.5):
    metals = {'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 
              'NA', 'K', 'CA', 'LI', 'RB', 'CS', 'MG', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 
              'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}

    # Get the script paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    lib_file_path = os.path.join(script_dir, lib_file)

    # Initialize PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)
    
    # Dictionaries and sets to store results
    metal_coordination = {}
    residues_to_extract = set()
    original_order = []  # Maintain the original residue order

    # Step 1: Analyze coordination and track residues
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_key = (chain.id, residue.get_id())
                original_order.append(residue_key)
                
                # Check for metal atoms
                for atom in residue:
                    if atom.element in metals:
                        residues_to_extract.add(residue_key)
                        metal_coord = atom.coord
                        key = f"{atom.element}_{chain.id}_{residue.get_id()}"
                        metal_coordination[key] = {
                            'metal_position': metal_coord,
                            'coordinating_residues': []
                        }

                        # Find nearby residues within cutoff
                        for chain2 in model:
                            for residue2 in chain2:
                                for atom2 in residue2:
                                    distance = np.linalg.norm(metal_coord - atom2.coord)
                                    if distance <= distance_cutoff:
                                        coord_info = {
                                            'chain': chain2.id,
                                            'residue_number': residue2.get_id()[1],
                                            'residue_name': residue2.get_resname(),
                                            'atom_name': atom2.get_name().strip(),
                                            'distance': round(distance, 2),
                                            'atom_coord': atom2.coord,
                                            'element': atom2.element
                                        }
                                        metal_coordination[key]['coordinating_residues'].append(coord_info)
                                        residues_to_extract.add((chain2.id, residue2.get_id()))

    # Step 2: Rebuild structure, explicitly handling hydrogens
    new_structure = PDB.Structure.Structure('metal_site')
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    new_chain = PDB.Chain.Chain('A')
    new_model.add(new_chain)

    # Add residues in original order, including all hydrogens
    for chain_id, res_id in original_order:
        if (chain_id, res_id) in residues_to_extract:
            original_residue = structure[0][chain_id][res_id]
            new_residue = PDB.Residue.Residue(original_residue.id, original_residue.resname, original_residue.segid)

            # Add all atoms, avoiding the skipping of duplicate names
            for atom in original_residue.get_atoms():
                new_atom = PDB.Atom.Atom(
                    atom.get_name(), atom.coord, atom.bfactor, atom.occupancy,
                    atom.altloc, atom.fullname, atom.serial_number, element=atom.element
                )
                new_residue.add(new_atom)
            new_chain.add(new_residue)

    # Step 3: Save the new structure
    ref_output = "REFQM.pdb"
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(ref_output)

    # Step 4: Extract specific atoms and write to fixed_charges.dat
    target_atoms = {'N', 'CA', 'C', 'O'}
    with open(ref_output, 'r') as pdb_file, open('fixed_charges.dat', 'w') as output_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if atom_name in target_atoms:
                    atom_number = int(line[6:11].strip())
                    residue_name = line[17:20].strip()
                    output_file.write(f"{atom_number} {atom_name} {residue_name}\n")

    # Step 5: Process MOL2 and library files
    mol2_atom_types = extract_atom_types_from_mol2(mol2_file)
    write_pdb_with_mol2_types(new_structure, output_pdb, mol2_atom_types)

    lib_types = read_aminofb15_lib(lib_file_path)
    update_atom_types_from_library(output_pdb, ref_output, lib_types)

    return metal_coordination

#Extract non-standard residues into part_QM files and create separate files for standard residues linked to metals.
def extract_non_standard_residues_from_ref(ref_pdb, output_pdb="part_QM.pdb", output_xyz="part_QM.xyz", standard_residues=None):
    metals = {'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 
              'NA', 'K', 'CA', 'LI', 'RB', 'CS', 'MG', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 
              'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}

    element_converter = {
        'ZN': 'Zn', 'FE': 'Fe', 'CO': 'Co', 'RU': 'Ru', 'IR': 'Ir', 'PT': 'Pt', 'AU': 'Au', 'AG': 'Ag', 
        'CU': 'Cu', 'MG': 'Mg', 'MN': 'Mn', 'NI': 'Ni', 'PD': 'Pd', 'CD': 'Cd', 'HG': 'Hg', 'BR': 'Br',
        'CL': 'Cl', 'AL': 'Al', 'GA': 'Ga', 'IN': 'In', 'SB': 'Sb', 'TL': 'Tl', 'PB': 'Pb', 'BI': 'Bi',
        'AS': 'As', 'SE': 'Se', 'SR': 'Sr', 'MO': 'Mo', 'TC': 'Tc', 'RE': 'Re', 'OS': 'Os', 'RH': 'Rh', 
        'NA': 'Na', 'BA': 'Ba', 'SI': 'Si'
    }

    if standard_residues is None:
        standard_residues = {
            "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE", 
            "HIP", "HYP", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }

    # Store metal coordinates
    metal_coords = []
    # Store residues and their atoms
    standard_res_atoms = {}  # For standard residues
    nonstandard_res_atoms = []  # For non-standard residues
    residue_name_mapping = {}  # Maps original residue identifiers to new names

    # First pass: collect metal coordinates
    with open(ref_pdb, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip()
                if element.upper() in metals:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    metal_coords.append(np.array([x, y, z]))

    # Second pass: identify and store residues
    residue_type_count = {}
    with open(ref_pdb, 'r') as f:
        current_res = None
        current_res_atoms = []
        
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resnum = int(line[22:26])
                res_key = (resname, chain, resnum)
                
                # Check if this is a new residue
                if res_key != current_res:
                    # Process previous residue if exists
                    if current_res is not None and current_res_atoms:
                        is_near_metal = False
                        for atom_line in current_res_atoms:
                            x = float(atom_line[30:38])
                            y = float(atom_line[38:46])
                            z = float(atom_line[46:54])
                            coord = np.array([x, y, z])
                            
                            for metal_coord in metal_coords:
                                if np.linalg.norm(metal_coord - coord) <= 2.5:
                                    is_near_metal = True
                                    break
                            if is_near_metal:
                                break
                        
                        if is_near_metal:
                            if current_res[0] in standard_residues:
                                # Generate new residue name (e.g., HI1, CY1)
                                base_name = current_res[0][:2].upper()
                                if base_name not in residue_type_count:
                                    residue_type_count[base_name] = 1
                                else:
                                    residue_type_count[base_name] += 1
                                new_name = f"{base_name}{residue_type_count[base_name]}"
                                
                                # Store the mapping
                                residue_name_mapping[current_res] = new_name
                                
                                # Update residue name in atom lines
                                updated_atoms = []
                                for atom_line in current_res_atoms:
                                    updated_line = (atom_line[:17] + 
                                                  f"{new_name:<3}" +
                                                  atom_line[20:])
                                    updated_atoms.append(updated_line)
                                standard_res_atoms[current_res] = updated_atoms
                            else:
                                nonstandard_res_atoms.extend(current_res_atoms)
                    
                    current_res = res_key
                    current_res_atoms = []
                
                current_res_atoms.append(line)
        
        # Process last residue
        if current_res is not None and current_res_atoms:
            is_near_metal = False
            for atom_line in current_res_atoms:
                x = float(atom_line[30:38])
                y = float(atom_line[38:46])
                z = float(atom_line[46:54])
                coord = np.array([x, y, z])
                
                for metal_coord in metal_coords:
                    if np.linalg.norm(metal_coord - coord) <= 2.5:
                        is_near_metal = True
                        break
            
            if is_near_metal:
                if current_res[0] in standard_residues:
                    base_name = current_res[0][:2].upper()
                    if base_name not in residue_type_count:
                        residue_type_count[base_name] = 1
                    else:
                        residue_type_count[base_name] += 1
                    new_name = f"{base_name}{residue_type_count[base_name]}"
                    residue_name_mapping[current_res] = new_name
                    
                    updated_atoms = []
                    for atom_line in current_res_atoms:
                        updated_line = (atom_line[:17] + 
                                      f"{new_name:<3}" +
                                      atom_line[20:])
                        updated_atoms.append(updated_line)
                    standard_res_atoms[current_res] = updated_atoms
                else:
                    nonstandard_res_atoms.extend(current_res_atoms)

    # Write files as before...
    if nonstandard_res_atoms:
        with open(output_pdb, 'w') as f:
            for line in nonstandard_res_atoms:
                f.write(line)
        
        xyz_atoms = []
        for line in nonstandard_res_atoms:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()
            if not element:
                alpha_chars = ''.join(c for c in atom_name if c.isalpha())
                element = alpha_chars[:2].upper() if len(alpha_chars) >= 2 else alpha_chars[0].upper()
            element = element_converter.get(element.upper(), element)
            xyz_atoms.append((element, x, y, z))
        
        with open(output_xyz, 'w') as f:
            f.write(f"{len(xyz_atoms)}\n")
            f.write(f"Generated from {ref_pdb} - Non-standard residues\n")
            for element, x, y, z in xyz_atoms:
                f.write(f"{element:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}\n")

    # Write standard residues to separate files
    for res_key, atoms in standard_res_atoms.items():
        new_name = residue_name_mapping[res_key]
        
        # Write PDB file
        pdb_filename = f"{new_name}.pdb"
        with open(pdb_filename, 'w') as f:
            for line in atoms:
                f.write(line)
        
        # Write XYZ file
        xyz_atoms = []
        for line in atoms:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()
            if not element:
                alpha_chars = ''.join(c for c in atom_name if c.isalpha())
                element = alpha_chars[:2].upper() if len(alpha_chars) >= 2 else alpha_chars[0].upper()
            element = element_converter.get(element.upper(), element)
            xyz_atoms.append((element, x, y, z))
        
        xyz_filename = f"{new_name}.xyz"
        with open(xyz_filename, 'w') as f:
            f.write(f"{len(xyz_atoms)}\n")
            f.write(f"Generated from {ref_pdb} - {new_name} {res_key[1]}:{res_key[2]}\n")
            for element, x, y, z in xyz_atoms:
                f.write(f"{element:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}\n")

    return residue_name_mapping

#Generate a PDB file containing all residues except those in part_QM.pdb,
#with updated residue names for metal-coordinating residues.    
def generate_nonstand_pdb(input_pdb, part_qm_pdb, output_pdb="nonstand.pdb"):
    try:
        # First get the mapping from HI*.pdb, CY*.pdb files to know which residues to rename
        rename_mapping = {}
        for filename in os.listdir():
            if (filename.endswith('.pdb') and 
                len(filename) == 7 and  # e.g., "HI1.pdb"
                filename[2].isdigit()):  
                new_name = filename[:3]  # e.g., "HI1"
                
                # Read the first ATOM/HETATM line to get residue info
                with open(filename, 'r') as f:
                    for line in f:
                        if line.startswith(("ATOM", "HETATM")):
                            resnum = int(line[22:26])
                            # Store mapping by resnum
                            rename_mapping[resnum] = new_name
                            break
        
        # Read residues to exclude from part_QM.pdb
        exclude_residues = set()
        with open(part_qm_pdb, 'r') as qm_file:
            for line in qm_file:
                if line.startswith(("ATOM", "HETATM")):
                    resnum = int(line[22:26])
                    exclude_residues.add(resnum)
        
        # Process input PDB and write output
        current_resnum = None
        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith("END"):
                    continue
                    
                if line.startswith(("ATOM", "HETATM")):
                    resnum = int(line[22:26])
                    
                    # Skip if residue should be excluded
                    if resnum in exclude_residues:
                        continue
                    
                    # Check if this residue needs to be renamed
                    if resnum in rename_mapping:
                        new_name = rename_mapping[resnum]
                        # Replace residue name while keeping the rest of the line intact
                        line = line[:17] + f"{new_name:<3}" + line[20:]
                    
                    outfile.write(line)
                
                elif not line.startswith(("ATOM", "HETATM")):
                    outfile.write(line)
            
            outfile.write("END\n")
            
    except Exception as e:
        print(f"Error generating nonstand.pdb: {e}")
        raise

#Extract charges and atom types for atoms in PDB files from MOL2 file.    
def extract_qm_charges(ref_pdb, mol2_file, residue_name_mapping=None):
    if residue_name_mapping is None:
        residue_name_mapping = {}
    
    try:
        # First, identify all PDB files to process
        pdb_files = []
         
        # Add part_QM.pdb first (for non-standard residues)
        if os.path.exists("part_QM.pdb"):
            pdb_files.append(("part_QM.pdb", "charge_qm.dat"))
        
        # Add all renamed standard residue PDB files
        for key, new_name in residue_name_mapping.items():
            pdb_filename = f"{new_name}.pdb"
            if os.path.exists(pdb_filename):
                charge_filename = f"charge_{new_name}.dat"
                pdb_files.append((pdb_filename, charge_filename))
        # Open charges_all.dat for writing the mapping
        with open("charges_all.dat", 'w') as charges_all_file:
            # First, handle part_QM.pdb for non-standard residues
            if os.path.exists("part_QM.pdb"):
                charges_all_file.write("QM.mol2 charge_qm.dat\n")
            
            # Write mappings for renamed standard residues
            for key, new_name in residue_name_mapping.items():
                pdb_filename = f"{new_name}.pdb"
                if os.path.exists(pdb_filename):
                    charges_all_file.write(f"{new_name}.mol2 charge_{new_name}.dat\n")
        # Process each PDB file
        for pdb_file, charge_output in pdb_files:
            # Read atom IDs and residue information from PDB
            qm_atom_ids = []
            atom_residue_mapping = {}
            with open(pdb_file, 'r') as pdb:
                for line in pdb:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            atom_id = int(line[6:11].strip())
                            resname = line[17:20].strip()
                            chain = line[21:22].strip()
                            resnum = int(line[22:26])
                            
                            # Check if this residue was renamed
                            original_key = (resname, chain, resnum)
                            new_resname = residue_name_mapping.get(original_key, resname)
                            
                            qm_atom_ids.append(atom_id)
                            atom_residue_mapping[atom_id] = new_resname
                        except ValueError as e:
                            print(f"Warning: Could not parse atom ID from PDB line: {line.strip()}")
                            continue
            
            # Read charges and atom types from MOL2 file
            charges_and_types = {}
            is_atom_section = False
            with open(mol2_file, 'r') as mol2:
                for line in mol2:
                    if "@<TRIPOS>ATOM" in line:
                        is_atom_section = True
                        continue
                    elif "@<TRIPOS>" in line and "ATOM" not in line:
                        is_atom_section = False
                        continue
                    
                    if is_atom_section and line.strip():
                        try:
                            parts = line.split()
                            atom_id = int(parts[0])
                            if atom_id in qm_atom_ids:
                                charge = float(parts[-1])
                                atom_type = parts[5]  # Extract atom type from column 6
                                charges_and_types[atom_id] = {
                                    'charge': charge, 
                                    'atom_type': atom_type
                                }
                        except (ValueError, IndexError) as e:
                            print(f"Warning: Could not parse MOL2 line: {line.strip()}")
                            continue
            
            # Write charges to output file
            with open(charge_output, 'w') as charge_file:
                for atom_id in qm_atom_ids:
                    if atom_id in charges_and_types:
                        data = charges_and_types[atom_id]
                        charge_file.write(f"{data['charge']:.6f} {data['atom_type']}\n")
                    else:
                        print(f"Warning: No charge found for atom ID {atom_id}")
                 
    except Exception as e:
        print(f"Error extracting charges: {e}")
        raise
    
#Generate easyPARM_residues.dat file containing information about renamed standard residues.
#Format: standard_residue.pdb standard_residue.mol2 standard_residue_name
#Only includes residues that were originally standard residues (like HID, CYS, etc.)
def generate_easyparm_residues(output_file="easyPARM_residues.dat"):
    # Define standard residues
    standard_residues = {
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
        "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    }

    try:
        # Find all special residue PDB files that came from standard residues
        special_residues = []
        for filename in os.listdir():
            if (filename.endswith('.pdb') and
                len(filename) == 7 and  # e.g., "HI1.pdb"
                filename[2].isdigit()):

                # Read the file to check if it was originally a standard residue
                with open(filename, 'r') as f:
                    original_residue = None
                    for line in f:
                        if line.startswith(("ATOM", "HETATM")):
                            # Get the base name of the original residue (first two letters)
                            base_name = line[17:19].strip()
                            # Check if any standard residue starts with these letters
                            is_standard = any(std_res.startswith(base_name)
                                            for std_res in standard_residues)
                            if is_standard:
                                # Get new residue name (e.g., "HI1")
                                new_name = filename[:3]

                                # Create the three columns
                                pdb_file = f"{new_name}.pdb"
                                mol2_file = f"{new_name}.mol2"
                                residue_name = new_name

                                special_residues.append((pdb_file, mol2_file, residue_name))
                            break

        # Sort the list for consistent output
        special_residues.sort()

        # Write to output file
        with open(output_file, 'w') as f:
            for pdb, mol2, name in special_residues:
                f.write(f"{pdb} {mol2} {name}\n")

    except Exception as e:
        print(f"Error generating {output_file}: {e}")
        raise

#Generate a PDB file containing all residues except those in part_QM.pdb
def generate_easyPARM_nonstand_pdb(input_pdb, part_qm_pdb, output_pdb="easynonstands.pdb"):
    try:
        # Read residues to exclude from part_QM.pdb
        exclude_residues = set()
        with open(part_qm_pdb, 'r') as qm_file:
            for line in qm_file:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resnum = line[22:26].strip()
                    exclude_residues.add((resname, chain, resnum))
        
        
        # Track the current residue being processed
        current_exclude = None
        
        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                # Skip END line
                if line.startswith("END"):
                    continue
                
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resnum = line[22:26].strip()
                    
                    # Check if this is a new residue
                    current_residue = (resname, chain, resnum)
                    
                    # If this residue is in exclude list, mark it to skip all atoms
                    if current_residue in exclude_residues:
                        current_exclude = current_residue
                        continue
                    
                    # If we were previously excluding and this is a new residue, reset
                    if current_exclude and current_residue != current_exclude:
                        current_exclude = None
                    
                    # Write the line if not currently excluding
                    if current_exclude is None:
                        outfile.write(line)
                
                # Always write non-residue lines except END
                elif not line.startswith(("ATOM", "HETATM")):
                    outfile.write(line)
        
    except Exception as e:
        print(f"Error generating nonstand.pdb: {e}")
        raise

def main():
    try:
        input_pdb = "metalloprotein_easyPARM.pdb"
        mol2_file = "NEW_COMPLEX.mol2"
        output_mol2 = "COMPLEX_updated.mol2"
        part_qm_pdb = "part_QM.pdb"
        part_qm_xyz = "part_QM.xyz"
        qm_pdb = "QM.pdb"
        lib_file = "amber_refernce.lib"
        nonstand_pdb = "nonstand.pdb"
        charge_file = "charge_qm.dat"
        easynonstand_pdb = "easynonstands.pdb"
        # Analyze metal site and generate files
        metal_coordination = analyze_and_extract_metal_site(input_pdb, mol2_file, lib_file, qm_pdb)

        # Update MOL2 file with QM.pdb atom types
        update_mol2_with_qm_types(mol2_file, qm_pdb, output_mol2)

        # Extract non-standard residues
        residue_name_mapping = extract_non_standard_residues_from_ref("REFQM.pdb", part_qm_pdb, part_qm_xyz)

        # Generate nonstand.pdb excluding the metal and ligand residues
        generate_nonstand_pdb(input_pdb, part_qm_pdb, nonstand_pdb)
 
        # Extract charges for QM region
        extract_qm_charges("REFQM.pdb", "NEW_COMPLEX.mol2", residue_name_mapping) 
        # After running extract_non_standard_residues_from_ref and generate_nonstand_pdb
        generate_easyparm_residues()
        generate_easyPARM_nonstand_pdb(input_pdb, part_qm_pdb, easynonstand_pdb)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"Error processing files: {e}")
        raise

if __name__ == "__main__":
    main()

