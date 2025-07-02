#!/usr/bin/env python3
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

import os
import numpy as np
import re
from collections import defaultdict
from collections import deque
import sys

class MolecularAnalyzer:
    # Define metal elements as a class attribute
    METAL_ELEMENTS = {'LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE',
                      'CO', 'NI', 'CU', 'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU',
                      'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM',
                      'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W',
                      'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'FR', 'RA', 'AC',
                      'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO',
                      'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', 'RG', 'CN', 'NH', 'FL', 'MC',
                      'LV', 'TS', 'OG'}
    def __init__(self):
        self.atoms = []  # List of atoms with their properties
        self.bonds = []  # List of bonds
        self.bond_dict = {}  # Dictionary for fast bond lookups
        self.distances = {}  # Dictionary to store atom-atom distances
        self.angles = {}  # Dictionary to store atom-atom-atom angles
        self.rings = []  # Store all detected rings
        self.aromatic_rings = [] # Store all detected aromatic rings
        self.ring_systems = {
            'isolated': [],  # Rings not connected to other rings
            'fused': [],     # Rings sharing at least two consecutive atoms (like naphthalene)
            'bridged': [],   # Rings connected by bridge atoms (like biphenyl)
            'spiro': []      # Rings connected by a single atom
        }
        
    #Read a mol2 file and extract atomic information and connectivity.
    def read_mol2(self, mol2_file):
        self.mol2_file = mol2_file  # Store the filename 
        with open(mol2_file, 'r') as f:
            lines = f.readlines()
        
        # Parse sections of mol2 file
        section = None
        atom_section = False
        bond_section = False
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            # Check for section headers
            if line.startswith('@<TRIPOS>'):
                section = line[9:].lower()
                atom_section = (section == 'atom')
                bond_section = (section == 'bond')
                continue
            
            # Parse atom section
            if atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5]
                    
                    # Extract element from atom_name
                    element = re.match(r'([A-Za-z]+)', atom_name).group(1)
                    
                    self.atoms.append({
                        'id': atom_id,
                        'name': atom_name,
                        'element': element,
                        'coords': [x, y, z],
                        'mol2_type': atom_type,
                        'connections': [],
                        'aromatic': False,
                        'in_rings': [],
                        'ring_junction': False,  # Is atom part of more than one ring
                        'bridge_atom': False,    # Is atom a bridge between rings
                        'spiro_atom': False,     # Is atom a spiro junction
                        'functional_group': None # Attached functional group
                    })
            
            # Parse bond section
            if bond_section:
                parts = line.split()
                if len(parts) >= 4:
                    bond_id = int(parts[0])
                    atom1_id = int(parts[1])
                    atom2_id = int(parts[2])
                    bond_type = parts[3]
                    bond = {
                        'id': bond_id,
                        'atom1': atom1_id,
                        'atom2': atom2_id,
                        'type': bond_type,
                        'original_type': bond_type,
                        'aromatic': bond_type == 'ar',
                        'in_rings': []  # Track which rings this bond is part of
                    }
                    self.bonds.append(bond)
                    
                    # Add bond to bond_dict for fast lookup
                    key = frozenset([atom1_id, atom2_id])
                    self.bond_dict[key] = bond
                    
                    # Add connection information to atoms
                    self.atoms[atom1_id-1]['connections'].append(atom2_id)
                    self.atoms[atom2_id-1]['connections'].append(atom1_id)
    
    #Read distance file with format: id1 id2 distance like atom1 atom2 distance_value
    def read_distance_file(self, distance_file):
        if not os.path.exists(distance_file):
            print(f"Warning: Distance file {distance_file} not found.")
            return
            
        with open(distance_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    distance = float(parts[2])
                    
                    # Store distances both ways for easy lookup
                    key = tuple(sorted([atom1_id, atom2_id]))
                    self.distances[key] = distance
    
    #Read angle file with format: id1 id2 id3 angle like atom1 atom2 atom3 angle_value
    def read_angle_file(self, angle_file):
        if not os.path.exists(angle_file):
            print(f"Warning: Angle file {angle_file} not found.")
            return
            
        with open(angle_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 4:
                    atom1_id = int(parts[0])
                    atom2_id = int(parts[1])
                    atom3_id = int(parts[2])
                    angle = float(parts[3])
                    
                    # Store the angle
                    key = (atom1_id, atom2_id, atom3_id)
                    self.angles[key] = angle
    
    #Determine hybridization state for each atom based on connectivity and geometry.
    def determine_hybridization(self):
        for i, atom in enumerate(self.atoms):
            atom_id = atom['id']
            element = atom['element'].upper()
            connections = atom['connections']
            conn_count = len(connections)
            
            # Default hybridization state
            hybridization = 'unknown' # This is for atom such as metals because they follow the ULS for atom type not GAFF
            
            # Basic hybridization rules
            if element == 'C':
                if conn_count == 6:
                    hybridization = 'sp3'
                elif conn_count == 5:
                    hybridization = 'sp3'
                elif conn_count == 4:
                    hybridization = 'sp3'
                elif conn_count == 3:
                    hybridization = 'sp2'
                elif conn_count == 2:
                    hybridization = 'sp'
                    
            elif element == 'N':
                if conn_count == 4:
                    hybridization = 'sp3'
                elif conn_count == 3:
                    hybridization = 'sp2'  # Assume sp2 for 3-connected N in rings 
                elif conn_count == 2:
                    hybridization = 'sp2'  # Common in heterocycles
                elif conn_count == 1:
                    hybridization = 'sp'
                    
            elif element == 'O':
                if conn_count == 2:
                    hybridization = 'sp3'
                    # Check angle if available
                    if len(connections) == 2:
                        angle_key1 = (connections[0], atom_id, connections[1])
                        angle_key2 = (connections[1], atom_id, connections[0])

                        if angle_key1 in self.angles and self.angles[angle_key1] > 110:
                            hybridization = 'sp2'
                        elif angle_key2 in self.angles and self.angles[angle_key2] > 110:
                            hybridization = 'sp2'
                elif conn_count == 3:
                    hybridization = 'sp3'
                elif conn_count == 1:
                    hybridization = 'sp2'
                    
            elif element in ['F', 'CL', 'BR', 'I']:
                hybridization = 'sp3'
            elif element == 'S':
                if conn_count == 2:
                    hybridization = 'sp3'
                    if len(connections) == 2:
                        angle_key1 = (connections[0], atom_id, connections[1])
                        angle_key2 = (connections[1], atom_id, connections[0])
                        
                        if angle_key1 in self.angles and self.angles[angle_key1] > 110:
                            hybridization = 'sp2'
                        elif angle_key2 in self.angles and self.angles[angle_key2] > 110:
                            hybridization = 'sp2'
                elif conn_count in [3, 4]:
                    hybridization = 'sp3'
                elif conn_count == 1:
                    hybridization = 'sp3'
                    
            elif element == 'P':
                if conn_count == 2:
                    hybridization = 'sp2'  # P with 3 connections typically sp3 with lone pair
                elif conn_count == 3:
                    hybridization = 'sp3'  # P with 3 connections typically sp3 with lone pair
                elif conn_count == 4:
                    hybridization = 'sp3'  # Tetrahedral arrangement
                elif conn_count == 5:
                    hybridization = 'sp3d'  # Trigonal bipyramidal
                elif conn_count == 6:
                    hybridization = 'sp3d2'  # Octahedral arrangement
                # If P shows unusual connectivity, at least assign a default
                else:
                    hybridization = 'sp3'  # Default for P if connection count is unusual

            elif element == 'B':
                if conn_count == 3:
                    hybridization = 'sp2'
                elif conn_count == 4:
                    hybridization = 'sp3'
                                
            elif element == 'H':
                hybridization = 's'
                
            # Store the determined hybridization
            self.atoms[i]['hybridization'] = hybridization

    #Ring detection with max_cycle_size limit and deduplication.
    def identify_rings(self):
        # Create adjacency list for the molecular graph
        adjacency = defaultdict(list)
        for bond in self.bonds:
            adjacency[bond['atom1']].append(bond['atom2'])
            adjacency[bond['atom2']].append(bond['atom1'])
        
        # Find all simple cycles using modified BFS approach
        all_rings = []
        visited_atom_pairs = set()
        
        def find_cycles(start_atom, max_cycle_size=12):
            # For each atom, try finding cycles through all its neighbors
            processed_neighbors = set()
            
            for neighbor in adjacency[start_atom]:
                # Avoid processing the same neighbor pair twice
                if neighbor in processed_neighbors:
                    continue
                
                processed_neighbors.add(neighbor)
                
                # Initialize BFS with the neighbor
                queue = deque([(neighbor, [start_atom, neighbor])])
                path_visited = set([(start_atom, neighbor)])
                
                while queue:
                    current, path = queue.popleft()
                    
                    # Don't proceed if path exceeds max size
                    if len(path) > max_cycle_size:
                        continue
                    
                    # Check all neighbors of current node
                    for next_atom in adjacency[current]:
                        # Found a cycle back to start
                        if next_atom == start_atom and len(path) > 2:
                            # Create a canonical representation of the ring
                            ring = path.copy()
                            ring_key = frozenset(ring)
                            
                            # Only add if we haven't seen this ring before
                            if ring_key not in visited_atom_pairs:
                                visited_atom_pairs.add(ring_key)
                                all_rings.append(ring)
                        
                        # Continue path if not creating a small cycle and not visited
                        elif (next_atom != start_atom and 
                              next_atom not in path and 
                              len(path) < max_cycle_size and
                              (current, next_atom) not in path_visited):
                            
                            path_visited.add((current, next_atom))
                            new_path = path + [next_atom]
                            queue.append((next_atom, new_path))
        
        # Find cycles starting from each atom
        for atom_id in range(1, len(self.atoms) + 1):
            find_cycles(atom_id)
        
        # Remove duplicate rings - although we've tried to avoid duplicates during generation
        unique_rings = []
        seen_rings = set()
        
        for ring in all_rings:
            # Create a canonical representation of the ring
            # Sort the ring atoms to ensure consistent representation
            canonical_ring = frozenset(ring)
            
            if canonical_ring not in seen_rings:
                seen_rings.add(canonical_ring)
                unique_rings.append(ring)
        
        # Sort rings by size for SSSR determination
        unique_rings.sort(key=len)
        
        # Implement Smallest Set of Smallest Rings (SSSR)
        essential_rings = []
        covered_bonds = set()
        
        for ring in unique_rings:
            # Convert ring to bond set
            ring_bonds = set()
            for i in range(len(ring)):
                bond = tuple(sorted([ring[i], ring[(i+1) % len(ring)]]))
                ring_bonds.add(bond)
            
            # Check if this ring covers new bonds
            if not ring_bonds.issubset(covered_bonds):
                essential_rings.append(ring)
                covered_bonds.update(ring_bonds)
        
        # Store the final ring set
        self.rings = essential_rings
        
        # Mark which rings each atom and bond belongs to
        for i, ring in enumerate(self.rings):
            for atom_id in ring:
                self.atoms[atom_id-1]['in_rings'].append(i)
            
            # Use bond_dict to mark bonds in this ring
            for j in range(len(ring)):
                atom1 = ring[j]
                atom2 = ring[(j+1) % len(ring)]
                key = frozenset([atom1, atom2])
                if key in self.bond_dict:
                    self.bond_dict[key]['in_rings'].append(i) 

    # DETECT AROMATICITY REGARDING CASES LIKE
    #1- Multiple heteroatoms in a single ring
    #2- Metal-coordinated rings and organometallic complexes
    #3- NHC (N-heterocyclic carbene) ligands
    #4- Proper application of Huckel's rule
    #5- Special handling for all types of heteroaromatic rings
    #6- Fused ring systems
    def detect_aromaticity(self):
        # Reset aromaticity flags in case mol2 has ar for aromatic
        for atom in self.atoms:
            atom['aromatic'] = False
        for bond in self.bonds:
            bond['aromatic'] = False
            if bond['type'] == 'ar':
                bond['type'] = '1'  # Reset aromatic bonds to single bonds initially
        
        # Step 1: Identify candidate rings (primarily 5 and 6-membered, but allow others too if it is available)
        candidate_rings = []
        for ring_idx, ring in enumerate(self.rings):
            # Most common aromatic rings are 5 or 6-membered, but we'll consider others
            if 3 <= len(ring) <= 7:
                candidate_rings.append(ring_idx)
        
        # Step 2: Create a list of metals and special handling elements
        # These need special treatment in aromatic systems detection
        metals = self.METAL_ELEMENTS
        # Group of atoms that can form aromatic ring
        group_15_16_atoms = ['N', 'P', 'AS', 'SB', 'O', 'S', 'SE', 'TE']
        
        # Step 3: Apply aromaticity criteria to each candidate ring
        aromatic_rings = []
        
        for ring_idx in candidate_rings:
            ring = self.rings[ring_idx]
            ring_size = len(ring)
            
            # Track heteroatoms that might participate in aromaticity despite appearing sp3 like presence of branches
            special_heteroatoms = set()
            
            # Track metal coordination sites within the ring
            metal_coordination_sites = set()
            
            # First pass: identify special atoms within the ring 
            for atom_id in ring:
                atom = self.atoms[atom_id-1]
                element = atom['element'].upper()
                
                # Skip hydrogen atoms. Consider only heavy atom
                if element == 'H':
                    continue
                
                # Check for connections to metals (for potential metal-coordinated aromatics)
                metal_connections = []
                for conn_id in atom['connections']:
                    conn_element = self.atoms[conn_id-1]['element'].upper()
                    if conn_element in metals:
                        metal_connections.append(conn_id)
                
                # If connected to metals, mark as coordination site
                if metal_connections:
                    metal_coordination_sites.add(atom_id)
                    special_heteroatoms.add(atom_id)  # Metal-coordinated atoms can still be part of aromatics
                
                # Special case for common heteroatoms (expanded list)
                if element in group_15_16_atoms + ['B']:
                    # Count connections excluding hydrogen
                    heavy_connections = sum(1 for conn in atom['connections'] 
                                        if self.atoms[conn-1]['element'].upper() != 'H')
                    
                    # Group 16 elements (S, O, Se, Te)
                    if element in ['O', 'S', 'SE', 'TE']:
                        if heavy_connections == 2:
                            special_heteroatoms.add(atom_id)
                        # Special case for sulfonyl or sulfoxide in aromatic rings
                        elif element == 'S' and heavy_connections > 2:
                            # Check if additional connections are exocyclic
                            ring_connections = sum(1 for conn in atom['connections'] 
                                               if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                            if ring_connections == 2:  # Two connections in ring, others outside
                                special_heteroatoms.add(atom_id)
                    
                    # Group 15 elements (N, P, As, Sb) - with special handling for metal coordination
                    elif element in ['N', 'P', 'AS', 'SB']:
                        # Always consider it capable of aromaticity if metal-coordinated
                        if metal_connections:
                            special_heteroatoms.add(atom_id)
                        elif heavy_connections <= 3:
                            # Check if additional connections (beyond 2) are exocyclic
                            if heavy_connections > 2:
                                ring_connections = sum(1 for conn in atom['connections'] 
                                                   if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                                if ring_connections >= 1:  # Allow even one connection in ring
                                    special_heteroatoms.add(atom_id)
                            else:
                                special_heteroatoms.add(atom_id)
                    
                    # Special case for boron (can be aromatic in certain systems)
                    elif element == 'B':
                        special_heteroatoms.add(atom_id)  # Always consider boron
            
            # Second pass: evaluate aromaticity with special handling for metal coordination
            # Count pi electrons and check if atoms can participate in aromaticity
            pi_electron_count = 0
            all_atoms_can_participate = True
            
            # Check if this ring has metal coordination - special case for NHC and similar ligands
            has_metal_coordination = len(metal_coordination_sites) > 0
            
            for atom_id in ring:
                atom = self.atoms[atom_id-1]
                element = atom['element'].upper()
                
                # Skip hydrogen atoms
                if element == 'H':
                    continue
                    
                # Check if atom has appropriate hybridization or special cases
                # Include metal-coordinated atoms and special heteroatoms
                if (atom['hybridization'] not in ['sp2', 'sp'] and 
                    atom_id not in special_heteroatoms and
                    atom_id not in metal_coordination_sites):
                    
                    # Special case: force check GAFF type if available (sometimes hybridization is misdetected)
                    if 'gaff_type' in atom and atom['gaff_type'] in ['n2', 'nd', 'nh', 'nb', 'nc', 'ne']:
                        # These GAFF types are typically planar and can participate in aromatics
                        pass
                    else:
                        all_atoms_can_participate = False
                        break
                
                # Count pi electrons - with special handling for metal coordination
                if element == 'C':
                    # Carbon contributes 1 pi electron unless it's a carbene (coordinated to metal)
                    if atom_id in metal_coordination_sites:
                        # Carbene carbon connected to metals in NHC - special case
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        if heavy_connections == 3 and any(self.atoms[conn-1]['element'].upper() in metals 
                                                       for conn in atom['connections']):
                            pi_electron_count += 0  # Neutral carbene carbon donating to metal
                        else:
                            pi_electron_count += 1
                    else:
                        pi_electron_count += 1
                
                elif element == 'N':
                    # For nitrogen, depends on connectivity, metal coordination, and position
                    if atom_id in metal_coordination_sites:
                        # For metal-coordinated N, only contribute 1 electron (the other forms the M-N bond)
                        pi_electron_count += 1
                    else:
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        
                        if heavy_connections == 2:
                            # N with 2 connections (pyrrole-like) contributes 2 pi electrons
                            pi_electron_count += 2
                        elif heavy_connections == 3:
                            # N with 3 connections - check configurations
                            ring_connections = sum(1 for conn in atom['connections'] 
                                                if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                            
                            if ring_connections == 2:  # N with exocyclic connection
                                pi_electron_count += 1  # Typically contributes 1 like in imidazole
                            else:
                                pi_electron_count += 1  # Fully in-ring with 3 connections
                
                elif element == 'O':
                    # Oxygen typically contributes 2 pi electrons
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 1  # If coordinated to metal
                    else:
                        pi_electron_count += 2
                
                elif element == 'S':
                    # Sulfur handling
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 1
                    else:
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        
                        if heavy_connections == 2:
                            pi_electron_count += 2
                        elif heavy_connections > 2:
                            # For sulfonyl/sulfoxide
                            ring_connections = sum(1 for conn in atom['connections'] 
                                               if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                            if ring_connections == 2:
                                pi_electron_count += 2
                
                elif element == 'P':
                    # Phosphorus handling
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 1
                    else:
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        
                        if heavy_connections == 2:
                            pi_electron_count += 2
                        elif heavy_connections == 3:
                            ring_connections = sum(1 for conn in atom['connections'] 
                                               if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                            if ring_connections == 2:
                                pi_electron_count += 2  # Phosphole-like
                            else:
                                pi_electron_count += 1
                        elif heavy_connections == 4:
                                pi_electron_count += 1  # Contributes 1 pi electron when double-bonded
                
                elif element in ['AS', 'SB']:
                    # Arsenic and antimony similar to phosphorus
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 1
                    else:
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        
                        if heavy_connections == 2:
                            pi_electron_count += 2
                        elif heavy_connections == 3:
                            ring_connections = sum(1 for conn in atom['connections'] 
                                               if conn in ring and self.atoms[conn-1]['element'].upper() != 'H')
                            if ring_connections == 2:
                                pi_electron_count += 2
                            else:
                                pi_electron_count += 1
                
                elif element in ['SE', 'TE']:
                    # Selenium and tellurium similar to sulfur
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 1
                    else:
                        pi_electron_count += 2
                
                elif element == 'B':
                    # Boron handling
                    if atom_id in metal_coordination_sites:
                        pi_electron_count += 0  # Depends on specific coordination
                    else:
                        heavy_connections = sum(1 for conn in atom['connections'] 
                                             if self.atoms[conn-1]['element'].upper() != 'H')
                        
                        if heavy_connections == 3:
                            pi_electron_count += 0  # Neutral boron often accepts electrons
                        elif heavy_connections == 2:
                            pi_electron_count += 1
            
            # Special case for metal-coordinated rings like NHC ligands
            # These should be considered aromatic even if they don't strictly follow Huckel's rule
            # Imidazole-type ligands with metal coordination are usually aromatic, hydrogen of N is replaced by the metal
            if has_metal_coordination and ring_size == 5 and all_atoms_can_participate:
                # Check if this is an imidazole-type structure (has 2 nitrogens)
                n_count = sum(1 for atom_id in ring 
                           if self.atoms[atom_id-1]['element'].upper() == 'N')
                
                if n_count >= 2:
                    # This is likely an NHC or similar ligand - mark as aromatic
                    aromatic_rings.append(ring_idx)
                    continue  # Skip the Huckel rule check
            
            # Apply Huckel's rule: (4n+2) pi electrons for aromatic systems
            if all_atoms_can_participate:
                # Calculate n based on pi_electron_count = (4n+2) --> n = (pi_electron_count - 2)/4
                n_value = (pi_electron_count - 2) / 4
                
                # Check if n is an integer or very close to an integer
                is_huckel_compliant = abs(n_value - round(n_value)) < 0.1 and n_value >= 0
                
                # Special handling for 5 and 6-membered rings
                if ring_size == 6:
                    # 6-membered rings typically need 6 pi electrons (n=1)
                    if pi_electron_count == 6 or (pi_electron_count == 10 and is_huckel_compliant):
                        aromatic_rings.append(ring_idx)
                
                elif ring_size == 5:
                    # 5-membered heteroaromatic rings typically have 6 pi electrons (n=1)
                    if (pi_electron_count == 6 or pi_electron_count == 10) and is_huckel_compliant:
                        aromatic_rings.append(ring_idx)
                    # Special case for coordinated systems with 5-membered rings
                    elif has_metal_coordination and 4 <= pi_electron_count <= 7:
                        aromatic_rings.append(ring_idx)
                
                # For other ring sizes, strictly apply Huckel's rule
                elif is_huckel_compliant:
                    aromatic_rings.append(ring_idx)
        
        # Step 4: Handle fused ring systems recursively
        original_count = -1
        new_count = len(aromatic_rings)
        
        # Keep expanding until no new rings are added
        while new_count > original_count:
            original_count = new_count
            
            additional_rings = []
            for ring_idx in aromatic_rings:
                ring = self.rings[ring_idx]
                
                for other_idx in candidate_rings:
                    # Skip rings already identified as aromatic
                    if other_idx in aromatic_rings or other_idx in additional_rings:
                        continue
                    
                    other_ring = self.rings[other_idx]
                    
                    # Check if rings are fused (share at least two atoms)
                    shared_atoms = set(ring) & set(other_ring)
                    if len(shared_atoms) >= 2:
                        # Check if shared atoms are already part of an aromatic system
                        shared_aromatic = all(self.atoms[atom_id-1]['hybridization'] in ['sp2', 'sp'] or
                                             any(self.atoms[conn-1]['element'].upper() in metals 
                                                for conn in self.atoms[atom_id-1]['connections'])
                                             for atom_id in shared_atoms
                                             if self.atoms[atom_id-1]['element'].upper() != 'H')
                        
                        if shared_aromatic:
                            # For fused systems, we need to be more flexible
                            additional_rings.append(other_idx)
            
            # Add new aromatic rings
            aromatic_rings.extend(additional_rings)
            new_count = len(aromatic_rings)
        
        # Step 5: Mark atoms in aromatic rings as aromatic
        for ring_idx in aromatic_rings:
            ring = self.rings[ring_idx]
            
            # Identify special cases within the ring
            special_heteroatoms = set()
            metal_coordination_sites = set()
            
            for atom_id in ring:
                atom = self.atoms[atom_id-1]
                element = atom['element'].upper()
                
                if element == 'H':
                    continue
                
                # Check for metal connections
                for conn_id in atom['connections']:
                    conn_element = self.atoms[conn_id-1]['element'].upper()
                    if conn_element in metals:
                        metal_coordination_sites.add(atom_id)
                
                # Identify special heteroatoms
                if element in group_15_16_atoms + ['B']:
                    heavy_connections = sum(1 for conn in atom['connections'] 
                                         if self.atoms[conn-1]['element'].upper() != 'H')
                    
                    # Consider all metal-coordinated atoms as special cases
                    if atom_id in metal_coordination_sites:
                        special_heteroatoms.add(atom_id)
                    
                    # Group 16 elements
                    elif element in ['O', 'S', 'SE', 'TE']:
                        if heavy_connections == 2:
                            special_heteroatoms.add(atom_id)
                        elif element == 'S' and heavy_connections > 2:
                            ring_connections = sum(1 for conn in atom['connections'] 
                                               if conn in ring and 
                                               self.atoms[conn-1]['element'].upper() != 'H')
                            if ring_connections == 2:
                                special_heteroatoms.add(atom_id)
                    
                    # Group 15 elements
                    elif element in ['N', 'P', 'AS', 'SB']:
                        if heavy_connections <= 3:
                            if heavy_connections > 2:
                                ring_connections = sum(1 for conn in atom['connections'] 
                                                   if conn in ring and 
                                                   self.atoms[conn-1]['element'].upper() != 'H')
                                if ring_connections >= 1:
                                    special_heteroatoms.add(atom_id)
                            else:
                                special_heteroatoms.add(atom_id)
                    
                    # Boron
                    elif element == 'B':
                        special_heteroatoms.add(atom_id)
            
            # Mark atoms as aromatic
            for atom_id in ring:
                atom = self.atoms[atom_id-1]
                # Mark as aromatic if meets criteria or is connected to a metal
                if atom['element'].upper() != 'H' and \
                   (atom['hybridization'] in ['sp2', 'sp'] or 
                    atom_id in special_heteroatoms or
                    atom_id in metal_coordination_sites or
                    # Check GAFF type as fallback
                    ('gaff_type' in atom and atom['gaff_type'] in ['n2', 'nh', 'na', 'nb', 'nc', 'ne'])):
                    atom['aromatic'] = True
        
        # Step 6: Mark bonds in aromatic rings as aromatic without metal coordination
        for bond in self.bonds:
            atom1_id = bond['atom1']
            atom2_id = bond['atom2']
            
            # Skip bonds involving metals directly
            if any(self.atoms[atom_id-1]['element'].upper() in metals 
                  for atom_id in [atom1_id, atom2_id]):
                continue
            
            # Both atoms must be aromatic
            if not (self.atoms[atom1_id-1]['aromatic'] and self.atoms[atom2_id-1]['aromatic']):
                continue
            
            # Check if they share an aromatic ring
            atom1_rings = set(self.atoms[atom1_id-1]['in_rings'])
            atom2_rings = set(self.atoms[atom2_id-1]['in_rings'])
            common_rings = atom1_rings & atom2_rings
            
            if common_rings and any(ring_idx in aromatic_rings for ring_idx in common_rings):
                bond['aromatic'] = True
                bond['type'] = 'ar'
        
        # Store aromatic rings for reference
        self.aromatic_rings = aromatic_rings
    
    #Calculate the center of mass for a ring. needed for certain cases like the bridge atoms.
    def calculate_ring_com(self, ring_idx):
        ring = self.rings[ring_idx]
        x_sum, y_sum, z_sum = 0, 0, 0
        total_mass = 0
        
        for atom_id in ring:
            atom = self.atoms[atom_id-1]
            # Skip hydrogen atoms in COM calculation
            if atom['element'].upper() == 'H':
                continue
                
            # Use atomic mass values for the elements
            # Simplified approach using approximate atomic masses
            element = atom['element'].upper()
            mass = self.get_atomic_mass(element)
            
            x, y, z = atom['coords']
            x_sum += x * mass
            y_sum += y * mass
            z_sum += z * mass
            total_mass += mass
        
        # Avoid division by zero
        if total_mass == 0:
            return None
            
        com = [x_sum / total_mass, y_sum / total_mass, z_sum / total_mass]
        return com

    #Get the atomic mass for an element.
    def get_atomic_mass(self, element):
        # Dictionary of common atomic masses
        atomic_masses = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
            'P': 30.974, 'S': 32.065, 'CL': 35.453, 'BR': 79.904, 'I': 126.904,
            'B': 10.811, 'SI': 28.085, 'SE': 78.971, 'TE': 127.60,
            'AS': 74.922, 'SB': 121.760, 'LI': 6.941, 'NA': 22.990,
            'K': 39.098, 'RB': 85.468, 'CS': 132.905, 'FR': 223.0,
            'BE': 9.012, 'MG': 24.305, 'CA': 40.078, 'SR': 87.62,
            'BA': 137.327, 'RA': 226.0, 'SC': 44.956, 'TI': 47.867,
            'V': 50.942, 'CR': 51.996, 'MN': 54.938, 'FE': 55.845,
            'CO': 58.933, 'NI': 58.693, 'CU': 63.546, 'ZN': 65.38,
            'Y': 88.906, 'ZR': 91.224, 'NB': 92.906, 'MO': 95.95,
            'TC': 98.0, 'RU': 101.07, 'RH': 102.906, 'PD': 106.42,
            'AG': 107.868, 'CD': 112.411, 'HF': 178.49, 'TA': 180.948,
            'W': 183.84, 'RE': 186.207, 'OS': 190.23, 'IR': 192.217,
            'PT': 195.084, 'AU': 196.967, 'HG': 200.592, 'LU': 174.967,
            'LR': 262.0, 'LA': 138.905, 'AC': 227.0
        }
        
        # Return mass or default to carbon if not found
        return atomic_masses.get(element, 12.0)

    #Calculate the distance between two centers of mass.
    def calculate_com_distance(self, com1, com2):
        return np.linalg.norm(np.array(com1) - np.array(com2))

    #Calculate the COM differences between fused rings.
    def calculate_fused_ring_com_differences(self):
        # Store COM differences
        self.fused_ring_com_diffs = {}
        
        # Calculate COMs for all rings
        ring_coms = {}
        for i in range(len(self.rings)):
            ring_coms[i] = self.calculate_ring_com(i)
        
        # Calculate differences for fused rings
        for i in range(len(self.rings)):
            for j in range(i+1, len(self.rings)):
                # Check if rings are in the fused category
                if i in self.ring_systems['fused'] and j in self.ring_systems['fused']:
                    # Check if they share atoms (direct connection)
                    shared_atoms = set(self.rings[i]) & set(self.rings[j])
                    if len(shared_atoms) >= 2:
                        # Check if they share consecutive atoms (fused rings)
                        consecutive = False
                        ring_i_list = self.rings[i]
                        ring_j_list = self.rings[j]
                        
                        # Find shared atoms that are adjacent in at least one ring
                        for atom1 in shared_atoms:
                            for atom2 in shared_atoms:
                                if atom1 != atom2:
                                    # Check adjacency in ring i
                                    idx1_i = ring_i_list.index(atom1)
                                    idx2_i = ring_i_list.index(atom2)
                                    adjacent_i = (abs(idx1_i - idx2_i) == 1 or 
                                                abs(idx1_i - idx2_i) == len(ring_i_list) - 1)
                                    
                                    # Check adjacency in ring j
                                    idx1_j = ring_j_list.index(atom1)
                                    idx2_j = ring_j_list.index(atom2)
                                    adjacent_j = (abs(idx1_j - idx2_j) == 1 or 
                                                abs(idx1_j - idx2_j) == len(ring_j_list) - 1)
                                    
                                    if adjacent_i or adjacent_j:
                                        consecutive = True
                                        break
                        
                        if consecutive and ring_coms[i] is not None and ring_coms[j] is not None:
                            # Calculate distance between COMs
                            distance = self.calculate_com_distance(ring_coms[i], ring_coms[j])
                            self.fused_ring_com_diffs[(i, j)] = {
                                'distance': distance,
                                'com1': ring_coms[i],
                                'com2': ring_coms[j]
                            }
                            
    #Classification of rings as isolated, fused, bridged, or spiro.
    def classify_ring_systems(self):
        # First, identify atoms that are part of multiple rings
        for atom in self.atoms:
            if len(atom['in_rings']) > 1:
                atom['ring_junction'] = True
        
        # Create a graph where nodes are rings and edges indicate connections
        ring_graph = defaultdict(set)
        
        # Identify connections between rings
        ring_connections = {}  # Store the atoms connecting each pair of rings
        
        # Find shared atoms between rings
        for i in range(len(self.rings)):
            for j in range(i+1, len(self.rings)):
                shared_atoms = set(self.rings[i]) & set(self.rings[j])
                if shared_atoms:
                    # Store connection
                    ring_graph[i].add(j)
                    ring_graph[j].add(i)
                    ring_connections[(i, j)] = list(shared_atoms)
                    ring_connections[(j, i)] = list(shared_atoms)
        
        # Find connected components (ring systems)
        visited_rings = set()
        ring_systems = []
        
        def find_ring_system(ring_idx):
            system = []
            queue = [ring_idx]
            visited_rings.add(ring_idx)
            
            while queue:
                current = queue.pop(0)
                system.append(current)
                
                for neighbor in ring_graph[current]:
                    if neighbor not in visited_rings:
                        visited_rings.add(neighbor)
                        queue.append(neighbor)
            
            return system
        
        # Reset ring system classifications
        self.ring_systems = {
            'isolated': [],  # Rings not connected to other rings
            'fused': [],     # Rings sharing at least two consecutive atoms
            'bridged': [],   # Rings connected by bridge atoms
            'spiro': []      # Rings connected by a single atom
        }
        
        # Find all ring systems
        for i in range(len(self.rings)):
            if i not in visited_rings:
                system = find_ring_system(i)
                if system:
                    ring_systems.append(system)
        
        # Classify each ring and ring system
        for system in ring_systems:
            if len(system) == 1:
                # Isolated ring
                self.ring_systems['isolated'].append(system[0])
            else:
                # System with multiple rings, classify connections
                system_connection_types = set()
                
                for i in range(len(system)):
                    for j in range(i+1, len(system)):
                        ring_i = system[i]
                        ring_j = system[j]
                        
                        if ring_j in ring_graph[ring_i]:
                            # These rings are connected, classify the connection
                            shared_atoms = ring_connections[(ring_i, ring_j)]
                            
                            if len(shared_atoms) == 1:
                                # Spiro junction
                                system_connection_types.add('spiro')
                                atom_id = shared_atoms[0]
                                self.atoms[atom_id-1]['spiro_atom'] = True
                            
                            elif len(shared_atoms) >= 2:
                                # Check if atoms are consecutive in both rings
                                consecutive = False
                                
                                # Get rings as ordered lists to check for adjacency
                                ring_i_list = self.rings[ring_i]
                                ring_j_list = self.rings[ring_j]
                                
                                # Find shared atoms that are adjacent in at least one ring
                                for atom1 in shared_atoms:
                                    for atom2 in shared_atoms:
                                        if atom1 != atom2:
                                            # Check adjacency in ring i
                                            idx1_i = ring_i_list.index(atom1)
                                            idx2_i = ring_i_list.index(atom2)
                                            adjacent_i = (abs(idx1_i - idx2_i) == 1 or 
                                                        abs(idx1_i - idx2_i) == len(ring_i_list) - 1)
                                            
                                            # Check adjacency in ring j
                                            idx1_j = ring_j_list.index(atom1)
                                            idx2_j = ring_j_list.index(atom2)
                                            adjacent_j = (abs(idx1_j - idx2_j) == 1 or 
                                                        abs(idx1_j - idx2_j) == len(ring_j_list) - 1)
                                            
                                            if adjacent_i or adjacent_j:
                                                consecutive = True
                                                break
                                
                                if consecutive:
                                    # Fused rings (sharing adjacent atoms)
                                    system_connection_types.add('fused')
                                else:
                                    # Bridged rings (sharing non-adjacent atoms)
                                    system_connection_types.add('bridged')
                                    for atom_id in shared_atoms:
                                        self.atoms[atom_id-1]['bridge_atom'] = True
                
                # Categorize the entire system based on the types of connections
                if 'fused' in system_connection_types:
                    self.ring_systems['fused'].extend(system)
                elif 'bridged' in system_connection_types:
                    self.ring_systems['bridged'].extend(system)
                elif 'spiro' in system_connection_types:
                    self.ring_systems['spiro'].extend(system)
        
        # Calculate COM differences for fused rings
        if self.ring_systems['fused']:
            self.calculate_fused_ring_com_differences()

    #Identifies ring junction atoms where rings of different sizes meet (like indole).
    #Identifies bridge atoms specifically for biphenyl-like connections between aromatic rings in presence of metal that lead to form a five ring with the bridge atoms.
    def identify_bridge_atoms(self):
        
        # First reset all bridge atom and ring junction flags
        for atom in self.atoms:
            atom['bridge_atom'] = False
            atom['ring_junction'] = False
        
        # Calculate COMs for all rings if not already calculated
        ring_coms = {}
        for i in range(len(self.rings)):
            ring_coms[i] = self.calculate_ring_com(i)
        
        # Identify atoms that are part of multiple rings (potential ring junctions)
        for atom in self.atoms:
            if len(atom['in_rings']) > 1:
                # Mark as ring junction
                atom['ring_junction'] = True
                
                # Evaluate ring properties
                ring_sizes = [len(self.rings[ring_idx]) for ring_idx in atom['in_rings']]
                
                # Special case for junctions between rings of different sizes
                if len(set(ring_sizes)) > 1:  # Different sized rings (e.g., 5-6 fusion)
                    atom['bridge_atom'] = True  # Mark as bridge for GAFF assigning process
                
                # For rings of the same size, check COM distance
                elif len(set(ring_sizes)) == 1:
                    # Get all pairs of rings this atom is part of
                    ring_pairs = []
                    for i, ring_i in enumerate(atom['in_rings']):
                        for j in range(i+1, len(atom['in_rings'])):
                            ring_j = atom['in_rings'][j]
                            ring_pairs.append((ring_i, ring_j))
                    
                    # Check COM difference for each pair
                    for ring_i, ring_j in ring_pairs:
                        # Make sure we have valid COMs for both rings
                        if ring_coms[ring_i] is not None and ring_coms[ring_j] is not None:
                            # Calculate COM distance
                            com_difference = self.calculate_com_distance(ring_coms[ring_i], ring_coms[ring_j])
                            
                            # Use COM difference threshold to identify bridge atoms
                            if com_difference >= 3.0:  # 3.0 Ã… threshold for significant spatial separation
                                atom['bridge_atom'] = True
                                # Optional: Store the COM difference for reference
                                if not hasattr(atom, 'com_differences'):
                                    atom['com_differences'] = {}
                                atom['com_differences'][(ring_i, ring_j)] = com_difference
                                 
        # Identify true bridge atoms between separate ring systems 
        for bond in self.bonds:
            # Skip bonds that are part of rings
            if bond['in_rings']:
                continue
            atom1_id = bond['atom1']
            atom2_id = bond['atom2']
            atom1 = self.atoms[atom1_id - 1]
            atom2 = self.atoms[atom2_id - 1]
            # Both atoms must be in rings, and their ring sets must differ
            if atom1['in_rings'] and atom2['in_rings']:
                ring_set1 = set(atom1['in_rings'])
                ring_set2 = set(atom2['in_rings'])
                # Check if the rings are aromatic
                rings1_aromatic = any(ring_idx in self.aromatic_rings for ring_idx in atom1['in_rings'])
                rings2_aromatic = any(ring_idx in self.aromatic_rings for ring_idx in atom2['in_rings'])
                # If they're in different rings and connected rings are aromatic, mark as bridge atoms
                if not ring_set1.intersection(ring_set2) and rings1_aromatic and rings2_aromatic:
                    # These are true bridge atoms (like in biphenyl)
                    if atom1['element'].upper() == 'C' and atom1['hybridization'] == 'sp2':
                        atom1['bridge_atom'] = True
                    if atom2['element'].upper() == 'C' and atom2['hybridization'] == 'sp2':
                        atom2['bridge_atom'] = True
          
        # Second filter: exclude non-carbon atoms and terminal atoms from bridge atoms
        for atom in self.atoms:
            if atom['element'].upper() != 'C':
                atom['bridge_atom'] = False
            # Terminal atoms (with only one connection) can never be bridge atoms
            if len(atom['connections']) <= 1:
                atom['bridge_atom'] = False
         
        # Third Filter: exclude bridge atoms connected to metal atoms
        metal_elements = self.METAL_ELEMENTS 
        
        excluded_by_metal = []
        for atom_idx, atom in enumerate(self.atoms, 1):
            if atom['bridge_atom']:
                # Check all connected atoms
                for connected_atom_id in atom['connections']:
                    connected_atom = self.atoms[connected_atom_id - 1]
                    if connected_atom['element'].upper() in metal_elements:
                        atom['bridge_atom'] = False
                        excluded_by_metal.append((atom_idx, atom, connected_atom_id, connected_atom['element']))
                        break
                        
        # Final Filter: exclude bridges where connected atoms have identical ring sets (Fused Ring)
        # Step 1: First identify all bridge bonds
        bridge_bonds = []
        for bond in self.bonds:
            atom1_id = bond['atom1']
            atom2_id = bond['atom2']
            atom1 = self.atoms[atom1_id - 1]
            atom2 = self.atoms[atom2_id - 1]
            
            if atom1['bridge_atom'] and atom2['bridge_atom']:
                bridge_bonds.append((atom1_id, atom2_id))

        # Step 2: Classify each bond as true bridge or same-ring bridge
        true_bridge_bonds = []
        same_ring_bridge_bonds = []

        for atom1_id, atom2_id in bridge_bonds:
            atom1 = self.atoms[atom1_id - 1]
            atom2 = self.atoms[atom2_id - 1]
            
            # If atoms share exactly the same rings, it's a same-ring bridge
            if set(atom1['in_rings']) == set(atom2['in_rings']):
                same_ring_bridge_bonds.append((atom1_id, atom2_id))
            else:
                true_bridge_bonds.append((atom1_id, atom2_id))

        # Step 3: Create set of atoms that are part of true bridge bonds
        true_bridge_atoms = set()
        for atom1_id, atom2_id in true_bridge_bonds:
            true_bridge_atoms.add(atom1_id)
            true_bridge_atoms.add(atom2_id)

        # Step 4: For each same-ring bridge bond, exclude atoms that aren't in true_bridge_atoms
        bridge_to_exclude = set()
        for atom1_id, atom2_id in same_ring_bridge_bonds:
            if atom1_id not in true_bridge_atoms:
                bridge_to_exclude.add(atom1_id)
            
            if atom2_id not in true_bridge_atoms:
                bridge_to_exclude.add(atom2_id)

        # Apply exclusion
        for atom_id in bridge_to_exclude:
            self.atoms[atom_id - 1]['bridge_atom'] = False

    #Specialized function to detect truely biphenyl bridge atoms (seperated rings)
    def detect_biphenyl_bridge_atoms(self):
        # First ensure we have aromaticity information
        if not hasattr(self, 'aromatic_rings'):
            if hasattr(self, 'detect_aromaticity'):
                self.detect_aromaticity()
        
        # Track which atoms are involved in biphenyl-like connections
        biphenyl_bridge_atoms = set()
        
        # Track rings with bridges already assigned
        rings_with_cp = set()
        
        # Get ring sizes
        ring_sizes = {}
        for i in range(len(self.rings)):
            ring_sizes[i] = len(self.rings[i])
        
        # Find bridge bonds that connect two separate 6-membered aromatic rings
        for bond in self.bonds:
            # Skip bonds that are part of rings
            if bond['in_rings']:
                continue
            
            atom1_id = bond['atom1']
            atom2_id = bond['atom2']
            atom1 = self.atoms[atom1_id - 1]
            atom2 = self.atoms[atom2_id - 1]
            
            # One atom must be carbon and aromatic, the other can be carbon or nitrogen and aromatic
            if (atom1['element'].upper() in ['C','N'] and atom2['element'].upper() in ['N','C'] and
                atom1.get('aromatic', False) and atom2.get('aromatic', False)):
                
                # Both atoms must be in rings but not the same rings
                if atom1['in_rings'] and atom2['in_rings']:
                    # Find 6-membered aromatic rings for each atom
                    ring1 = None
                    ring2 = None
                    
                    for r in atom1['in_rings']:
                        if r in self.aromatic_rings and ring_sizes.get(r) == 6:
                            ring1 = r
                            break
                    
                    for r in atom2['in_rings']:
                        if r in self.aromatic_rings and ring_sizes.get(r) == 6:
                            ring2 = r
                            break
                    
                    # Skip if either atom is not in a 6-membered aromatic ring
                    if ring1 is None or ring2 is None or ring1 == ring2:
                        continue
                    
                    # Count heavy atom connections for each
                    heavy_conn1 = sum(1 for conn in atom1['connections'] 
                                   if self.atoms[conn-1]['element'].upper() != 'H')
                    heavy_conn2 = sum(1 for conn in atom2['connections'] 
                                   if self.atoms[conn-1]['element'].upper() != 'H')
                    
                    # For biphenyl bridge atoms, typically exactly 3 heavy atom connections
                    if heavy_conn1 == 3 and heavy_conn2 == 3:
                        # For each atom in the bridge, assign appropriate GAFF type
                        # Determine bridge type first (for the bond as a whole)
                        if ring1 in rings_with_cp or ring2 in rings_with_cp:
                            bridge_type = 'cq'  # Subsequent bridges get 'cq'
                        else:
                            bridge_type = 'cp'  # First bridge gets 'cp'
                            # Mark these rings as having a bridge
                            rings_with_cp.add(ring1)
                            rings_with_cp.add(ring2)

                        # Now assign types to individual atoms, overriding for nitrogen
                        for atom_id, atom in [(atom1_id, atom1), (atom2_id, atom2)]:
                            if atom['element'].upper() == 'N':
                                # Nitrogen in aromatic bridge always gets 'nb'
                                self.atoms[atom_id - 1]['gaff_type'] = 'nb'
                            else:
                                # Carbon atoms follow the bridge type
                                self.atoms[atom_id - 1]['gaff_type'] = bridge_type
                            
                            self.atoms[atom_id - 1]['bridge_atom'] = True
                            biphenyl_bridge_atoms.add(atom_id)                        
                        # Mark these rings as having a bridge
                        rings_with_cp.add(ring1)
                        rings_with_cp.add(ring2)
        
        return biphenyl_bridge_atoms
 
    #Identify common functional groups in the molecule.
    def identify_functional_groups(self):
        for atom in self.atoms:
            # Skip atoms in rings as we're looking for functional groups attached to rings
            if atom['in_rings']:
                continue
                
            atom_id = atom['id']
            element = atom['element'].upper()
            
            # Check for common functional groups based on patterns
            if element == 'O':
                # Check for carbonyl, alcohol, ether
                if len(atom['connections']) == 1:
                    conn_atom = self.atoms[atom['connections'][0]-1]
                    if conn_atom['element'].upper() == 'C':
                        # Check if carbon is double-bonded to this oxygen using bond_dict
                        key = frozenset([atom_id, conn_atom['id']])
                        if key in self.bond_dict and self.bond_dict[key]['type'] in ['2', 'ar']:
                            atom['functional_group'] = 'carbonyl'
                        else:
                            atom['functional_group'] = 'alcohol/ether'
            
            elif element == 'N':
                # Check for amine, amide, nitro
                if len(atom['connections']) <= 3:
                    conn_atoms = [self.atoms[conn_id-1] for conn_id in atom['connections']]
                    
                    # Check for amide (N connected to carbonyl C)
                    for conn_atom in conn_atoms:
                        if conn_atom['element'].upper() == 'C':
                            # Check if carbon is connected to oxygen with double bond
                            for conn_id2 in conn_atom['connections']:
                                if conn_id2 != atom_id and self.atoms[conn_id2-1]['element'].upper() == 'O':
                                    key = frozenset([conn_atom['id'], conn_id2])
                                    if key in self.bond_dict and self.bond_dict[key]['type'] == '2':
                                        atom['functional_group'] = 'amide'
                                        break
                    
                    # If not amide, it's likely an amine
                    if not atom['functional_group']:
                        atom['functional_group'] = 'amine'

    #Propagate conjugated types (e/f suffix) through a conjugated system.
    #Works for C, N, P and other elements in conjugated systems.
    #Handles both sp and sp2 atoms in conjugation.
    def propagate_conjugated_types(self, start_atom_id, start_type_suffix):
        # Initialize tracking variables if they don't exist
        if not hasattr(self, 'conjugated_types'):
            self.conjugated_types = {}
        
        # Store the full type prefix (element-dependent) for each atom
        element_prefixes = {
            'C': 'c',  # Carbon -> ce/cf or c1/cg
            'N': 'n',  # Nitrogen -> ne/nf
            'P': 'p',  # Phosphorus -> pe/pf
        }
        
        visited = set([start_atom_id])
        queue = [(start_atom_id, start_type_suffix)]
        
        # Assign the starting atom's type
        start_atom = self.atoms[start_atom_id-1]
        start_element = start_atom['element'].upper()
        start_hybridization = start_atom['hybridization']
        
        if start_element in element_prefixes:
            prefix = element_prefixes[start_element]
            # Handle sp hybridized atoms specially
            if start_hybridization == 'sp':
                # sp atoms in conjugation are 'cg' for carbon
                self.conjugated_types[start_atom_id] = 'cg' if prefix == 'c' else prefix + '1'
            else:
                self.conjugated_types[start_atom_id] = prefix + start_type_suffix
        
        while queue:
            current_id, current_suffix = queue.pop(0)
            current_atom = self.atoms[current_id-1]
            current_hybridization = current_atom['hybridization']
            connections = current_atom['connections']
            
            for conn_id in connections:
                if conn_id in visited:
                    continue
                    
                conn_atom = self.atoms[conn_id-1]
                conn_element = conn_atom['element'].upper()
                conn_hybridization = conn_atom['hybridization']
                
                # Only process atoms that could be part of conjugation and aren't aromatic
                if conn_element in element_prefixes and not conn_atom['aromatic']:
                    key = frozenset([current_id, conn_id])
                    if key not in self.bond_dict:
                        continue
                        
                    bond_type = self.bond_dict[key]['type']
                    
                    # Handle sp hybridized atoms
                    if conn_hybridization == 'sp':
                        prefix = element_prefixes[conn_element]
                        # For carbon in conjugation, use 'cg'
                        self.conjugated_types[conn_id] = 'cg' if prefix == 'c' else prefix + '1'
                        visited.add(conn_id)
                        
                        # Continue propagation from this sp atom
                        queue.append((conn_id, current_suffix))  # Use same suffix for propagation
                    
                    # Handle sp2 hybridized atoms
                    elif conn_hybridization == 'sp2':
                        prefix = element_prefixes[conn_element]
                        
                        # Determine suffix based on bond type
                        if bond_type == '1':  # Single bond
                            next_suffix = current_suffix  # Same suffix for single bonds
                        else:  # Double bond
                            next_suffix = 'f' if current_suffix == 'e' else 'e'  # Opposite suffix
                        
                        self.conjugated_types[conn_id] = prefix + next_suffix
                        visited.add(conn_id)
                        queue.append((conn_id, next_suffix))
    
    #Determine the chemical environment of each atom.
    def determine_chemical_environment(self):
        for i, atom in enumerate(self.atoms):
            element = atom['element'].upper()
            hybridization = atom['hybridization']
            is_aromatic = atom['aromatic']
            is_bridge = atom['bridge_atom']
            connections = atom['connections']
            
            # Initialize environment description
            environment = []
            
            # Basic atom type
            env_type = f"{element} ({hybridization})"
            if is_aromatic:
                env_type += " - aromatic"
                
            # Changed: Use "bridge carbon" for carbon atoms in multiple rings instead of "ring junction"
            if element == 'C' and atom['ring_junction'] and len(atom['in_rings']) >= 2:
                env_type += " - bridge carbon"
            elif atom['ring_junction']:  # Keep "ring junction" for non-carbon atoms
                env_type += " - ring junction"
                
            if atom['spiro_atom']:
                env_type += " - spiro junction"
            
            environment.append(env_type)
            
            # Determine connected atoms and their types
            connected_atoms_info = []
            for conn_id in connections:
                conn_atom = self.atoms[conn_id-1]
                conn_element = conn_atom['element'].upper()
                conn_aromatic = conn_atom['aromatic']
                conn_bridge = conn_atom['bridge_atom']
                
                conn_type = conn_element
                if conn_aromatic:
                    conn_type += " (aromatic)"
                if conn_bridge:
                    conn_type += " (bridge)"
                
                connected_atoms_info.append(conn_type)
            
            # Count element types
            element_counts = defaultdict(int)
            for conn_type in connected_atoms_info:
                element_counts[conn_type] += 1
            
            # Build connections string
            conn_str = []
            for elem_type, count in element_counts.items():
                if count == 1:
                    conn_str.append(f"{elem_type}")
                else:
                    conn_str.append(f"{count}-{elem_type}")
            
            if conn_str:
                environment.append("connected to " + ", ".join(conn_str))

            atom['conn_str'] = conn_str
            
            # Add ring membership with improved accuracy
            if atom['in_rings']:
                rings_info = []
                for ring_idx in atom['in_rings']:
                    ring = self.rings[ring_idx]
                    ring_type = "unknown"
                    
                    # Classify ring
                    if ring_idx in self.ring_systems['isolated']:
                        ring_type = "isolated"
                    elif ring_idx in self.ring_systems['fused']:
                        ring_type = "fused"
                    elif ring_idx in self.ring_systems['bridged']:
                        ring_type = "bridged"
                    elif ring_idx in self.ring_systems['spiro']:
                        ring_type = "spiro"
                    
                    # Check if ring is aromatic
                    ring_aromatic = all(self.atoms[atom_id-1]['aromatic'] for atom_id in ring 
                                       if self.atoms[atom_id-1]['element'].upper() != 'H')
                    
                    if ring_aromatic:
                        ring_type += " aromatic"
                    
                    rings_info.append(f"{len(ring)}-membered {ring_type} ring")
                
                environment.append("part of " + ", ".join(rings_info))
            
            # Add functional group info if any
            if atom['functional_group']:
                environment.append(f"functional group: {atom['functional_group']}")
            
            # Store the final environment string
            self.atoms[i]['environment'] = " - ".join(environment)


    #Assigns GAFF atom types to each atom based on its chemical properties:
    #element, hybridization, aromaticity, chemical environment, connectivity, etc.
    def assign_gaff_atom_types(self):
        # Track aromatic ring atoms and heterocyclic rings for cc/ca differentiation
        biphenyl_bridge_atoms = self.detect_biphenyl_bridge_atoms()
        aromatic_ring_atoms = set()
        heterocyclic_rings = set()
        
        # Track ring sizes and heterocycle status
        for ring_idx, ring in enumerate(self.rings):
            # Check if ring is heterocyclic (contains non-carbon atoms)
            non_carbon_atoms = [atom_id for atom_id in ring 
                                if self.atoms[atom_id-1]['element'].upper() != 'C' and 
                                   self.atoms[atom_id-1]['element'].upper() != 'H']
            
            if non_carbon_atoms:
                heterocyclic_rings.add(ring_idx) 
                # Add all aromatic atoms from this ring
                for atom_id in ring:
                    if self.atoms[atom_id-1]['aromatic']:
                        aromatic_ring_atoms.add(atom_id)
        
        for i, atom in enumerate(self.atoms):
            element = atom['element'].upper()
            hybridization = atom['hybridization']
            is_aromatic = atom['aromatic']
            connections = atom['connections']
            is_bridge = atom['bridge_atom']
            is_in_rings = atom['in_rings']
            atom_id = atom['id']
            if atom_id in biphenyl_bridge_atoms:
            # These were already set to 'cp' and/or 'cq' in detect_biphenyl_bridge_atoms
                continue
            
            # Initialize with None (for metals and non-standard elements)
            gaff_type = None
            
            # Handle hydrogen atoms - depends on what they're connected to
            if element == 'H':
                if len(connections) == 1:
                    attached_atom = self.atoms[connections[0] - 1]
                    attached_element = attached_atom['element'].upper()
                    if attached_element == 'C':
                        attached_hybridization = attached_atom['hybridization']
                        # Count EWGs (N, O, F, Cl, Br, I) attached to the carbon
                        ewg_count = sum(1 for conn_id in attached_atom['connections']
                                       if self.atoms[conn_id - 1]['element'].upper() in ['N', 'O', 'F', 'CL', 'BR', 'I', 'S', 'P'])
                        
                        if attached_hybridization == 'sp3':
                            # Aliphatic sp3 carbon
                            if ewg_count == 1:
                                gaff_type = 'h1'
                            elif ewg_count == 2:
                                gaff_type = 'h2'
                            elif ewg_count == 3:
                                gaff_type = 'h3'
                            else:
                                gaff_type = 'hc'  # Default for spÃ‚Â³ carbon
                        else:
                            # Non-sp3 carbon (sp2 or sp)
                            if attached_atom['aromatic']:
                                if ewg_count == 0:
                                    gaff_type = 'ha'
                                elif ewg_count == 1:
                                    gaff_type = 'h4'
                                elif ewg_count == 2:
                                    gaff_type = 'h5'
                                else:
                                    gaff_type = 'ha'  # Default for aromatic
                            else:
                                # Non-aromatic, non-sp3 carbon
                                if ewg_count == 0:
                                    gaff_type = 'hc'
                                elif ewg_count == 1:
                                    gaff_type = 'h4'
                                elif ewg_count == 2:
                                    gaff_type = 'h5'
                                else:
                                    gaff_type = 'hc'  # Default for non-aromatic
                    elif attached_element == 'N':
                        gaff_type = 'hn'
                    elif attached_element == 'O':
                        gaff_type = 'ho'
                    elif attached_element == 'P':
                        gaff_type = 'hp'
                    elif attached_element == 'S':
                        gaff_type = 'hs'
                    else:
                            gaff_type = 'ha'  # Generic hydrogen 
            # Handle carbon atoms
            elif element == 'C':
                if hybridization == 'sp3':
                    # Check for ring membership and size
                    if is_in_rings:
                        ring_sizes = [len(self.rings[ring_idx]) for ring_idx in is_in_rings]
                        min_ring_size = min(ring_sizes)
                        # Check if any atom in the ring is a metal
                        has_metal_in_ring = False
                        for ring_idx in is_in_rings:
                            ring = self.rings[ring_idx]
                            for ring_atom_id in ring:
                                ring_atom = self.atoms[ring_atom_id-1]
                                ring_element = ring_atom['element'].upper()
                                # Define metal elements 
                                metal_elements = self.METAL_ELEMENTS 
                                if ring_element in metal_elements:
                                    has_metal_in_ring = True
                                    break
                            if has_metal_in_ring:
                                break
                        
                        if has_metal_in_ring:
                            gaff_type = 'c3'
                        else:
                          
                            if min_ring_size == 3:
                                gaff_type = 'cx'  # sp3 C in three-membered ring
                            elif min_ring_size == 4:
                                gaff_type = 'cy'  # sp3 C in four-membered ring
                            elif min_ring_size == 5:
                                gaff_type = 'c5'  # sp3 C in five-membered ring
                            elif min_ring_size == 6:
                                gaff_type = 'c6'  # sp3 C in six-membered ring
                            else:
                                gaff_type = 'c3'  # Generic sp3 C in larger rings
                    else:
                        gaff_type = 'c3'  # Generic sp3 C
                        
                elif hybridization == 'sp2':
                    # Check for carbonyl carbon
                    has_double_bond_to_o = False
                    has_double_bond_to_s = False
                    
                    for conn_id in connections:
                        conn_atom = self.atoms[conn_id-1]
                        key = frozenset([atom_id, conn_id])
                        if key in self.bond_dict:
                            bond_type = self.bond_dict[key]['type']
                            if conn_atom['element'].upper() == 'O' and bond_type == '2':
                                has_double_bond_to_o = True
                            elif conn_atom['element'].upper() == 'S' and bond_type == '2':
                                has_double_bond_to_s = True
                    
                    if has_double_bond_to_o:
                        gaff_type = 'c'  # sp2 C in C=O
                    elif has_double_bond_to_s:
                        gaff_type = 'cs'  # sp2 C in C=S
                    elif is_aromatic:
                        # Distinction between bridge atoms and ring junctions
                        if atom['bridge_atom'] and len(is_in_rings) == 2:
                                #Either true bridge atoms (like in biphenyl) or junctions between rings of different sizes
                            gaff_type = 'cp'  # sp2 C bridging aromatic rings or at 5-6 ring junctions
                        elif atom['ring_junction']:
                            # Ring junction atoms of same-sized rings (like in naphthalene)
                            gaff_type = 'ca'  # Inner sp2 C at ring junctions in same-sized rings
                        elif any(len(self.rings[ring_idx]) == 5 for ring_idx in is_in_rings):

                            has_heteroatom_connection = False
                            # Get the elements of connected atoms directly
                            element_list = []
                            for conn_id in connections:
                                connected_atom = self.atoms[conn_id-1]
                                element_list.append(connected_atom['element'])

                            # Filter out hydrogen atoms for determining atom types
                            non_h_elements = [elem for elem in element_list if elem != 'H']

                            # Define when to use cd vs cc based on specific element combinations
                            if len(non_h_elements) >= 2:
                                # cd is only when it's between two C, or C-N, or C-P
                                if all(elem == 'C' for elem in non_h_elements) or \
                                   (non_h_elements.count('C') == 1 and non_h_elements.count('N') == 1) or \
                                   (non_h_elements.count('C') == 1 and non_h_elements.count('P') == 1):
                                    has_heteroatom_connection = False  # Use cd
                                else:
                                    has_heteroatom_connection = True   # Use cc
                            gaff_type = 'cc' if has_heteroatom_connection else 'cd'

                        else:
                            gaff_type = 'ca'  # sp2 C in six-membered aromatic rings
                    
                    elif is_in_rings:
                        ring_size = min(len(self.rings[ring_idx]) for ring_idx in is_in_rings)
                        # Check if any atom in the ring is a metal
                        has_metal_in_ring = False
                        for ring_idx in is_in_rings:
                            ring = self.rings[ring_idx]
                            for ring_atom_id in ring:
                                ring_atom = self.atoms[ring_atom_id-1]
                                ring_element = ring_atom['element'].upper()
                                # Define metal elements 
                                metal_elements = self.METAL_ELEMENTS 
                                if ring_element in metal_elements:
                                    has_metal_in_ring = True
                                    break
                            if has_metal_in_ring:
                                break
                        
                        if has_metal_in_ring:
                            gaff_type = 'c'
                        else:
                            if ring_size == 3:
                                gaff_type = 'cu'  # sp2 C in three-membered ring
                            elif ring_size == 4:
                                gaff_type = 'cv'  # sp2 C in four-membered ring
                            elif ring_size == 5:
                                # Get all directly connected elements from connection list
                                has_heteroatom_connection = False

                                # Get the elements of connected atoms directly
                                element_list = []
                                for conn_id in connections:
                                    connected_atom = self.atoms[conn_id-1]
                                    element_list.append(connected_atom['element'])

                                non_h_elements = [elem for elem in element_list if elem != 'H']

                                if len(non_h_elements) >= 2:
                                    # cd is only when it's between two C, or C-N, or C-P
                                    if all(elem == 'C' for elem in non_h_elements) or \
                                       (non_h_elements.count('C') == 1 and non_h_elements.count('N') == 1) or \
                                       (non_h_elements.count('C') == 1 and non_h_elements.count('P') == 1):
                                        has_heteroatom_connection = False  # Use cd
                                    else:
                                        has_heteroatom_connection = True   # Use cc

                                # If connected to any heteroatom, use cc, otherwise cd
                                gaff_type = 'cc' if has_heteroatom_connection else 'cd'
                            elif ring_size >= 6:
                                    # Get all directly connected elements from connection list
                                    has_heteroatom_connection = False

                                    # Get the elements of connected atoms directly
                                    element_list = []
                                    for conn_id in connections:
                                        connected_atom = self.atoms[conn_id-1]
                                        element_list.append(connected_atom['element'])

                                    non_h_elements = [elem for elem in element_list if elem != 'H']

                                    if len(non_h_elements) >= 2:
                                        # cd is only when it's between two C, or C-N, or C-P
                                        if (non_h_elements.count('C') == 1 and non_h_elements.count('N') == 1) or \
                                           (non_h_elements.count('C') == 1 and non_h_elements.count('P') == 1):
                                            has_heteroatom_connection = True  # Use cd
                                        else:
                                            has_heteroatom_connection = False   # Use cc

                                    # If connected to any heteroatom, use cd, otherwise cc
                                    gaff_type = 'cd' if has_heteroatom_connection else 'cc'
                    else:
                        # Check if this is part of a guanidine group
                        connected_to_n = sum(1 for conn_id in connections
                                          if self.atoms[conn_id-1]['element'].upper() == 'N')
                        if connected_to_n >= 2:
                            # This might be a guanidine group
                            has_double_bond_to_n = False
                            for conn_id in connections:
                                if self.atoms[conn_id-1]['element'].upper() == 'N':
                                    key = frozenset([atom_id, conn_id])
                                    if key in self.bond_dict and self.bond_dict[key]['type'] == '2':
                                        has_double_bond_to_n = True
                                        break
                            if has_double_bond_to_n and connected_to_n >= 2:
                                gaff_type = 'cz'  # sp2 carbon in guanidine group
                            else:
                                gaff_type = 'c2'  # Generic sp2 C
                        else:
                            # Check if carbon is part of a conjugated system
                            is_in_conjugated_system = False
                            
                            # We'll need to know which atoms are already assigned
                            # This could be stored in a class variable or dictionary
                            if not hasattr(self, 'conjugated_types'):
                                self.conjugated_types = {}
                            
                            # First, identify if this atom is part of a conjugated system
                            if not is_aromatic:  # Skip aromatic atoms
                                for conn_id in connections:
                                    conn_atom = self.atoms[conn_id-1]
                                    if conn_atom['element'].upper() in ['C', 'N', 'P'] and conn_atom['hybridization'] == 'sp2':
                                        key = frozenset([atom_id, conn_id])
                                        if key in self.bond_dict and self.bond_dict[key]['type'] == '1':
                                            # Found a single bond to another sp2 carbon - potential conjugation
                                            is_in_conjugated_system = True
                                            break

                            # For a carbon atom detected to be in a conjugated system:
                            if is_in_conjugated_system:
                                # If this atom is already part of an identified conjugated system
                                if atom_id in self.conjugated_types:
                                    gaff_type = self.conjugated_types[atom_id]
                                else:
                                    # Start a new conjugated system analysis
                                    # Choose a starting suffix (arbitrarily 'e')
                                    self.propagate_conjugated_types(atom_id, 'e')
                                    gaff_type = self.conjugated_types[atom_id]

                            else:
                                gaff_type = 'c2'  # Generic sp2 C

                elif hybridization == 'sp':
                    # Check if this sp carbon is part of a conjugated system
                    # A carbon is considered part of a conjugated system if:
                    # 1. It's connected to an sp2 carbon (which would be ce or cf)
                    # 2. It's also connected to another sp (cg) or sp2 (ce,cf) carbon

                    connected_to_sp2 = False
                    connected_to_sp_or_sp2 = 0

                    for conn_id in connections:
                        conn_atom = self.atoms[conn_id-1]
                        if conn_atom['element'].upper() == 'C':
                            if conn_atom['hybridization'] == 'sp2':
                                connected_to_sp2 = True
                                connected_to_sp_or_sp2 += 1
                            elif conn_atom['hybridization'] == 'sp':
                                connected_to_sp_or_sp2 += 1

                    # If this sp carbon is in a conjugated system (connected to both sp2 and another sp/sp2)
                    if connected_to_sp2 and connected_to_sp_or_sp2 >= 2:
                        gaff_type = 'cg'  # Inner sp C in conjugated systems
                    else:
                        gaff_type = 'c1'  # Generic sp C
                else:
                    gaff_type = 'c2'  # Generic carbon
 
            # Handle nitrogen atoms
            elif element == 'N':
                if hybridization == 'sp3':
                    # Count connections and hydrogens
                    h_count = sum(1 for conn_id in connections 
                                if self.atoms[conn_id-1]['element'].upper() == 'H')
                    
                    if len(connections) == 4:
                        if h_count == 0:
                            gaff_type = 'n4'  # Quaternary N+
                        elif h_count == 1:
                            gaff_type = 'nx'  # Tertiary N+ with one H
                        elif h_count == 2:
                            gaff_type = 'ny'  # Secondary N+ with two H
                        elif h_count == 3:
                            gaff_type = 'nz'  # Primary N+ with three H
                        else:  # h_count == 4
                            gaff_type = 'n+'  # NH4+
                    elif len(connections) == 3:
                        if h_count == 0:
                            gaff_type = 'n3'  # Tertiary amine N
                        elif h_count == 1:
                            gaff_type = 'n7'  # Secondary amine N with one H
                        elif h_count == 2:
                            gaff_type = 'n8'  # Primary amine N with two H
                        else:  # h_count == 3
                            gaff_type = 'n9'  # NH3
                    else:
                        # Check if it's amine connected to aromatic ring
                        connected_to_aromatic = any(self.atoms[conn_id-1]['aromatic'] for conn_id in connections)
                        if connected_to_aromatic:
                            if h_count == 0:
                                gaff_type = 'nh'  # Tertiary amine N connected to aromatic rings
                            elif h_count == 1:
                                gaff_type = 'nu'  # Secondary amine N with one H, connected to aromatic ring
                            else:  # h_count >= 2
                                gaff_type = 'nv'  # Primary amine N with two H, connected to aromatic ring
                        else:
                            gaff_type = 'n3'  # Default sp3 N
                
                elif hybridization == 'sp2':
                    if is_aromatic:
                        ring_size = min(len(self.rings[ring_idx]) for ring_idx in is_in_rings)
                        if ring_size == 5 :
                            if len(connections) == 2:
                                # Get all directly connected elements from connection list
                                has_heteroatom_connection = False
                                # Get the elements of connected atoms directly
                                element_list = []
                                for conn_id in connections:
                                    connected_atom = self.atoms[conn_id-1]
                                    element_list.append(connected_atom['element'])
                                # Check if any connected atom is not C or H
                                has_heteroatom_connection = any(elem not in ('C','P', 'H') for elem in element_list)
                                # If connected to any heteroatom, use cc, otherwise cd
                                gaff_type = 'nc' if has_heteroatom_connection else 'nd'
                            elif len(connections) == 3:
                                # Check if all bonds are of type '1' (single bonds) or 'ar' (aromatic bonds)
                                if len(is_in_rings) == 1:
                                    gaff_type = 'na'
                                elif len(is_in_rings) in [2,3]:

                                    all_single_bonds = True
                                    for conn_id in connections:
                                        key = frozenset([atom_id, conn_id])
                                        if key in self.bond_dict:
                                            # Use original_type from the bond dictionary
                                            original_bond_type = self.bond_dict[key].get('original_type', self.bond_dict[key]['type'])
                                            if original_bond_type != '1':
                                                all_single_bonds = False
                                                break
                                
                                    # Assign 'na' only if all bonds were originally single bonds ('1')
                                    if all_single_bonds:
                                        gaff_type = 'na'
                                    else:
                                        gaff_type = 'nb'  # Default for N with non-single bonds                            
                        else:
                            gaff_type = 'nb'  # sp2 N for pure aromatic ring systems
                    else:
                        # Check for nitro group
                        connected_to_o = sum(1 for conn_id in connections 
                                           if self.atoms[conn_id-1]['element'].upper() == 'O')
                        if connected_to_o >= 2:
                            gaff_type = 'no'  # N in nitro group
                        elif len(connections) == 3:
                            connected_to_aromatic = any(self.atoms[conn_id-1]['aromatic'] for conn_id in connections)
                            if connected_to_aromatic:
                                gaff_type = 'na'  # sp2 N with three connected atoms
                            elif is_in_rings:
                                gaff_type = 'na'
                            else:
                                gaff_type = 'n3'

                        elif len(connections) == 2:
                            # Check if this is in amide group
                            ring_size = min(len(self.rings[ring_idx]) for ring_idx in is_in_rings)
                            if ring_size == 5:
                                # Get all directly connected elements from connection list
                                has_heteroatom_connection = False
                                # Get the elements of connected atoms directly
                                element_list = []
                                for conn_id in connections:
                                    connected_atom = self.atoms[conn_id-1]
                                    element_list.append(connected_atom['element'])
                                # Check if any connected atom is not C or H
                                has_heteroatom_connection = any(elem not in ('C','P', 'H') for elem in element_list)
                                # If connected to any heteroatom, use cc, otherwise cd
                                gaff_type = 'nc' if has_heteroatom_connection else 'nd'
                            else:
                             
                                h_count = sum(1 for conn_id in connections
                                            if self.atoms[conn_id-1]['element'].upper() == 'H')
                                
                                # Check for amide
                                connected_to_carbonyl = False
                                for conn_id in connections:
                                    conn_atom = self.atoms[conn_id-1]
                                    if conn_atom['element'].upper() == 'C':
                                        # Check if C is connected to O with double bond
                                        for c_conn_id in conn_atom['connections']:
                                            if self.atoms[c_conn_id-1]['element'].upper() == 'O':
                                                key = frozenset([conn_id, c_conn_id])
                                                if key in self.bond_dict and self.bond_dict[key]['type'] == '2':
                                                    connected_to_carbonyl = True
                                                    break
                                
                                # Check if N is part of a conjugated system (similar to C conjugation check)
                                is_in_conjugated_system = False
                                if not hasattr(self, 'conjugated_types'):
                                    self.conjugated_types = {}
                                
                                # Check if nitrogen is sp2 hybridized and part of a conjugated system
                                if hybridization == 'sp2' and not is_aromatic:
                                    for conn_id in connections:
                                        conn_atom = self.atoms[conn_id-1]
                                        # Look for connections to sp2 C or N atoms 
                                        if (conn_atom['element'].upper() in ['N', 'C'] and 
                                            conn_atom['hybridization'] == 'sp2'):
                                            key = frozenset([atom_id, conn_id])
                                            # Check if there's a single or double bond connection
                                            if key in self.bond_dict:
                                                is_in_conjugated_system = True
                                                break
                                
                                # Assign atom type based on findings
                                if is_in_conjugated_system:
                                    # If this atom is already part of an identified conjugated system
                                    if atom_id in self.conjugated_types:
                                        gaff_type = self.conjugated_types[atom_id]
                                    else:
                                        # Start a new conjugated system analysis with suffix 'e'
                                        self.propagate_conjugated_types(atom_id, 'e')
                                        gaff_type = self.conjugated_types[atom_id]
                                elif connected_to_carbonyl:
                                    if h_count == 0:
                                        gaff_type = 'n'    # Amide N with no H
                                    elif h_count == 1:
                                        gaff_type = 'ns'   # Amide N with one H
                                    else:  # h_count >= 2
                                        gaff_type = 'nt'   # Amide N with two H
                                else:
                                    # Default sp2 nitrogen with 2 connections that isn't in a conjugated system
                                    gaff_type = 'n2'   
 
                elif hybridization == 'sp':
                    gaff_type = 'n1'  # sp N (triple bonded)
                else:
                    gaff_type = 'n'

             
            elif element == 'O':
                if hybridization == 'sp3':
                    # Check for hydroxyl or ether/ester
                    is_hydroxyl = False
                    for conn_id in connections:
                        # Check connected atoms for hydrogens
                        conn_atom = self.atoms[conn_id-1]
                        if conn_atom['element'].upper() == 'H':
                            is_hydroxyl = True
                            break
                    
                    if is_hydroxyl:
                        gaff_type = 'oh'  # sp3 O in hydroxyl group
                    else:
                        # Check if this is water
                        if len(connections) == 2 and all(self.atoms[conn_id-1]['element'].upper() == 'H' for conn_id in connections):
                            gaff_type = 'ow'  # Oxygen in water
                        else:
                            gaff_type = 'os'  # sp3 O in ether and ester
                
                elif hybridization == 'sp2':
                    # Check for carbonyl or carboxylate
                    is_carbonyl = False
                    for conn_id in connections:
                        conn_atom = self.atoms[conn_id-1]
                        if conn_atom['element'].upper() == 'C':
                            # Check for double bond
                            key = frozenset([atom_id, conn_id])
                            if key in self.bond_dict and self.bond_dict[key]['type'] == '2':
                                is_carbonyl = True
                                break
                    
                    if is_carbonyl or len(connections) == 1:
                        gaff_type = 'o'  # sp2 O in C=O, COO-
                    else:
                        gaff_type = 'o'  # Default to ether/ester oxygen
                else:
                    gaff_type = 'os'
            
            # Handle other common elements
            elif element == 'S':
                # Count connections to determine if hypervalent sulfur
                if len(connections) >= 3:
                    if len(connections) == 3:
                        # Check for linked atom (sp2)
                        sp2_carbon_count = 0
                        for conn_id in connections:
                            connected_atom = self.atoms[conn_id-1]
                            if connected_atom['hybridization'] == 'sp2' and connected_atom['element'] in ['C', 'O','S', 'P']:
                                sp2_carbon_count += 1
                        
                        if sp2_carbon_count >= 2:
                            # Connected to one sp2 carbon, one sp2 oxygen  and one sp3 carbon
                            gaff_type = 'sx'  
                        else:
                            # Connected to one sp2 carbon and two sp3 carbons or all are sp3
                            gaff_type = 's4'  
                    else:  # 4 or more connections
                        # Count sp2 carbon connections
                        sp2_carbon_count = 0
                        for conn_id in connections:
                            connected_atom = self.atoms[conn_id-1]
                            if connected_atom['hybridization'] == 'sp2' and connected_atom['element'] == 'C':
                                sp2_carbon_count += 1
                        
                        if sp2_carbon_count >= 1:
                            # Connected to three sp2 atoms (including the vinyl group) and one sp3
                            gaff_type = 'sy'  # S with two double bonds to oxygen and connection to vinyl/sp2 carbon
                        else:
                            # Connected to two sp2 (the oxygens) and two sp3 (cyclopropyl groups)
                            gaff_type = 's6'  # S with two double bonds to oxygen and two single bonds to sp3 carbons                
                else:
                    if len(connections) == 1:
                        gaff_type = 's'  # Terminal S
                    elif len(connections) == 2:
                        # Check for double bond
                        has_double_bond = False
                        for conn_id in connections:
                            key = frozenset([atom_id, conn_id])
                            if key in self.bond_dict and self.bond_dict[key]['type'] == '2':
                                has_double_bond = True
                                break
                        
                        if has_double_bond:
                            gaff_type = 's2'  # S with at least one double bond
                        else:
                            # Check for thiol
                            has_h = any(self.atoms[conn_id-1]['element'].upper() == 'H' for conn_id in connections)
                            if has_h:
                                gaff_type = 'sh'  # S in thiol group
                            else:
                                gaff_type = 'ss'  # S in other sulfides
            
            elif element == 'P':
                if len(connections) >= 2 and not is_in_rings:
                    if len(connections) == 2:
                        # Check if P is part of a conjugated system
                        is_in_conjugated_system = False
                        
                        # We'll need to know which atoms are already assigned
                        # This could be stored in a class variable or dictionary
                        if not hasattr(self, 'conjugated_types'):
                            self.conjugated_types = {}
                        
                        # First, identify if this atom is part of a conjugated system
                        if not is_aromatic:  # Skip aromatic atoms
                            for conn_id in connections:
                                conn_atom = self.atoms[conn_id-1]
                                if conn_atom['element'].upper() in ['C', 'N', 'P'] and conn_atom['hybridization'] == 'sp2':
                                    key = frozenset([atom_id, conn_id])
                                    if key in self.bond_dict and self.bond_dict[key]['type'] == '1':
                                        # Found a single bond to another sp2 carbon - potential conjugation
                                        is_in_conjugated_system = True
                                        break

                        # For a P atom detected to be in a conjugated system:
                        if is_in_conjugated_system:
                            # If this atom is already part of an identified conjugated system
                            if atom_id in self.conjugated_types:
                                gaff_type = self.conjugated_types[atom_id]
                            else:
                                # Start a new conjugated system analysis
                                # Choose a starting suffix (arbitrarily 'e') (pe/pf atom type)
                                self.propagate_conjugated_types(atom_id, 'e')
                                gaff_type = self.conjugated_types[atom_id]
                        else:
                            gaff_type = 'p2'  # Phosphate with two connected atoms
                    elif len(connections) == 3:
                        # Check for hydrogen connections first
                        hydrogen_count = 0
                        sp2_carbon_count = 0
                        
                        for conn_id in connections:
                            connected_atom = self.atoms[conn_id-1]
                            if connected_atom['element'] == 'H':
                                hydrogen_count += 1
                            elif connected_atom['hybridization'] == 'sp2' and connected_atom['element'] == 'C':
                                sp2_carbon_count += 1
                        
                        if hydrogen_count > 0:
                            # P connected to at least one hydrogen atom
                            gaff_type = 'p3'
                        elif sp2_carbon_count >= 1:
                            # Connected to at least one sp2 carbon
                            gaff_type = 'px'
                        else:
                            # Connected to hydroxyl groups and sp3 carbon
                            gaff_type = 'p4'                        
                    else:  # 4 or more connections
                        has_double_bond_to_oxygen = False
                        connected_atoms = []

                        for conn_id in connections:
                            connected_atom = self.atoms[conn_id-1]
                            connected_atoms.append(connected_atom)

                            key = tuple(sorted([atom['id'], conn_id]))
                            if key in self.bond_dict and self.bond_dict[key]['type'] == '2' and connected_atom['element'] in ['O', 'S', 'C']:
                                has_double_bond_to_carbon = True

                        if has_double_bond_to_oxygen:
                            gaff_type = 'py'  # P with one double bonds to oxygen and connection to sp2 carbon
                        else:
                            gaff_type = 'p5'  # Phosphate with four connected atoms
                
                # Aromatic phosphorus 
                elif is_in_rings:
                    if is_aromatic:
                        gaff_type = 'pb'  # Sp2 P in pure aromatic systems
                    else:
                        # Checking ring atoms and not aromatic ring
                        has_heteroatom_connection = False

                        # Get the elements of connected atoms directly
                        element_list = []
                        for conn_id in connections:
                            connected_atom = self.atoms[conn_id-1]
                            element_list.append(connected_atom['element'])

                        # Check if any connected atom is not C or H
                        has_heteroatom_connection = any(elem not in ('C', 'H') for elem in element_list)

                        # If connected to any heteroatom, use pc, otherwise pd
                        gaff_type = 'pc' if has_heteroatom_connection else 'pd'
            
            # Handle halogens
            elif element == 'F':
                gaff_type = 'f'
            elif element == 'CL':
                gaff_type = 'cl'
            elif element == 'BR':
                gaff_type = 'br'
            elif element == 'I':
                gaff_type = 'i'
            elif element == 'B':
                gaff_type = 'b'
            
            # For metals and other elements not covered by GAFF, use the element symbol
            else:
                gaff_type = element.capitalize()
            
            # Store the determined GAFF type
            self.atoms[i]['gaff_type'] = gaff_type


        # Post-processing to handle double bonds between two 'cd' atoms
        cd_atoms = []
        for i, atom in enumerate(self.atoms):
            if atom.get('gaff_type') == 'cd':
                cd_atoms.append(atom['id'])
        # Check all pairs of 'cd' atoms to see if they have a double bond between them
        for i in range(len(cd_atoms)):
            for j in range(i+1, len(cd_atoms)):
                atom1_id = cd_atoms[i]
                atom2_id = cd_atoms[j]
                # Check if these atoms have a bond between them
                key = frozenset([atom1_id, atom2_id])
                if key in self.bond_dict:
                    # Get the bond type
                    #bond_type = self.bond_dict[key]['type']
                    original_bond_type = self.bond_dict[key].get('original_type', self.bond_dict[key]['type'])
                    # If it's a double bond (can be written as ar or 2), change one to 'cc'
                    if original_bond_type == '2':
                        # Check if either atom is connected to an 'nd' atom
                        atom1_connected_to_nd = False
                        atom2_connected_to_nd = False
                        
                        # Check atom1's connections
                        for bond_key in self.bond_dict:
                            if atom1_id in bond_key:
                                other_atom_id = next(id for id in bond_key if id != atom1_id)
                                if self.atoms[other_atom_id-1].get('gaff_type') == 'nd':
                                    atom1_connected_to_nd = True
                                    break
                        
                        # Check atom2's connections
                        for bond_key in self.bond_dict:
                            if atom2_id in bond_key:
                                other_atom_id = next(id for id in bond_key if id != atom2_id)
                                if self.atoms[other_atom_id-1].get('gaff_type') == 'nd':
                                    atom2_connected_to_nd = True
                                    break
                        
                        # Change the atom that is NOT connected to 'nd' to 'cc'
                        if atom1_connected_to_nd and not atom2_connected_to_nd:
                            self.atoms[atom2_id-1]['gaff_type'] = 'cc'
                        elif atom2_connected_to_nd and not atom1_connected_to_nd:
                            self.atoms[atom1_id-1]['gaff_type'] = 'cc'
                        # If both or neither are connected to 'nd', fall back to ID comparison
                        elif atom1_id < atom2_id:
                            self.atoms[atom1_id-1]['gaff_type'] = 'cc'
                        else:
                            self.atoms[atom2_id-1]['gaff_type'] = 'cc'
    #Perform full analysis of the molecular structure and assign the gaff atom type.
    def analyze(self, mol2_file, distance_file=None, angle_file=None):
        print(f"Analyzing {mol2_file}...")
        
        # Read input files
        self.read_mol2(mol2_file)
        if distance_file:
            self.read_distance_file(distance_file)
        if angle_file:
            self.read_angle_file(angle_file)
        
        # Identify rings and ring systems
        self.identify_rings()
        self.classify_ring_systems()
        self.identify_bridge_atoms()
        
        # Determine hybridization
        self.determine_hybridization()
        
        # Detect aromaticity
        self.detect_aromaticity()
        
        # Identify functional groups
        self.identify_functional_groups()
        
        # Determine chemical environments
        self.determine_chemical_environment()
        
        # Assign GAFF atom types
        self.assign_gaff_atom_types()
        
        # Print results
        self.print_results()
    
    #Get statistics about rings in the molecule.
    def get_ring_info(self):
        if not self.rings:
            return "No rings found"

        result = []
        result.append(f"Total rings: {len(self.rings)}")

        # Count by size
        ring_sizes = defaultdict(int)
        for ring in self.rings:
            ring_sizes[len(ring)] += 1

        result.append("Ring sizes:")
        for size, count in sorted(ring_sizes.items()):
            result.append(f"  {size}-membered rings: {count}")

        return "\n".join(result)

    #Print the analysis results about connectivity, hybridization, element and chemical Environment
    def print_results(self):
        print("\n=== MOLECULAR STRUCTURE ANALYSIS ===\n")
        print(f"Total atoms: {len(self.atoms)}")
        print(f"Total bonds: {len(self.bonds)}")
        
        # Print ring statistics
        print("\n--- RING STATISTICS ---")
        print(self.get_ring_info())
        
        print("\n--- ATOM PROPERTIES ---")
        with open("easyPARM_atomtype.dat", "w") as output: 
            for atom in self.atoms:
                print(f"\nAtom {atom['id']} ({atom['name']}):")
                print(f"  Element: {atom['element'].upper()}")
                print(f"  Hybridization: {atom['hybridization']}")
                print(f"  Aromatic: {'Yes' if atom['aromatic'] else 'No'}")
                print(f"  GAFF Type: {atom['gaff_type']}")
                output.writelines(f"{atom['gaff_type']}\n")
                print(f"  Chemical Environment: {atom['environment']}")
                print(f"  Connected to atoms: {', '.join(map(str, atom['connections']))}")
                if atom['in_rings']:
                    print(f"  Part of rings: {', '.join(str(r+1) for r in atom['in_rings'])}")

    #Generate a new MOL2 file with GAFF atom types replacing the original atom types.
    def generate_gaff_mol2(self, output_mol2="easyPARM.mol2"):
        
        try:
            # Process MOL2 file
            updated_mol2_lines = []
            with open(self.mol2_file, 'r') as input_file:
                is_atom_section = False
                atom_index = 0
                for line in input_file:
                    # Detect and modify the atom section
                    if line.startswith("@<TRIPOS>ATOM"):
                        is_atom_section = True
                        updated_mol2_lines.append(line)
                        atom_index = 0
                        continue
                    elif line.startswith("@<TRIPOS>"):
                        is_atom_section = False
                        updated_mol2_lines.append(line)
                        continue

                    if is_atom_section and line.strip():
                        # Match MOL2 atom line format using regex
                        atom_match = re.match(
                            r"(\s*\d+\s+)(\S+\s+)(-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+\s+)(\d+\s+\S+\s+)(-?\d+\.\d+)",
                            line
                        )
                        if atom_match and atom_index < len(self.atoms):
                            atom_id, atom_name, coords, atom_type, post_type, charge = atom_match.groups()
                            # Replace original atom type with GAFF atom type
                            gaff_type = self.atoms[atom_index]['gaff_type']
                            updated_line = f"{atom_id}{atom_name}{coords}{gaff_type:<7s}{post_type}{charge}\n"
                            updated_mol2_lines.append(updated_line)
                            atom_index += 1
                        else:
                            updated_mol2_lines.append(line)
                    else:
                        updated_mol2_lines.append(line)

            # Write the updated MOL2 content to the output file
            with open(output_mol2, 'w') as output_file:
                output_file.writelines(updated_mol2_lines)

            print(f"Successfully created GAFF mol2 file: {output_mol2}")
            return output_mol2

        except Exception as e:
            print(f"Error generating GAFF mol2 file: {e}")
            return None

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("Usage: python script.py <mol2_file> [distance_file] [angle_file] ")
        sys.exit(1)
    
    mol2_file = sys.argv[1]
    distance_file = sys.argv[2] if len(sys.argv) > 2 else None
    angle_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    analyzer = MolecularAnalyzer()
    analyzer.analyze(mol2_file, distance_file, angle_file)
    analyzer.generate_gaff_mol2()

