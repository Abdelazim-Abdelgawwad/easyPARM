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

from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set
import math
import os
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class FrcmodBond:
    atom1: str
    atom2: str
    force: float
    length: float

@dataclass
class FrcmodAngle:
    atom1: str
    atom2: str
    atom3: str
    force: float
    angle: float

@dataclass
class FrcmodDihedral:
    atom1: str
    atom2: str
    atom3: str
    atom4: str
    periodicity: int
    force: float
    phase: float
    divider: float

@dataclass
class FrcmodNonbond:
    atom_type: str
    radius: float
    well_depth: float

@dataclass
class CharmmBond:
    atom1: str
    atom2: str
    force: float
    length: float
    comment: str = ""

@dataclass
class CharmmAngle:
    atom1: str
    atom2: str
    atom3: str
    force: float
    angle: float
    comment: str = ""

@dataclass
class CharmmDihedral:
    atom1: str
    atom2: str
    atom3: str
    atom4: str
    force: float
    periodicity: int
    phase: float
    comment: str = ""

@dataclass
class CharmmNonbond:
    atom_type: str
    ignored: float  # Often 0.0 in CHARMM
    epsilon: float
    rmin_half: float
    comment: str = ""


#Converts molecular force field parameters between AMBER frcmod and CHARMM formats. First Guess    
class ParameterConverter:
    
    def __init__(self):
        # AMBER to CHARMM atom type mapping
        self.amber_to_charmm = {
            'c':   'CG2O1',
            'c1':  'CG1T1',
            'c2':  'CG2D1',
            'c3':  'CG331',
            'ca':  'CG2R61',
            'cp':  'CG2R61',
            'cq':  'CG2R67',
            'cc':  'CG2R57',
            'cd':  'CG2R51',
            'ce':  'CG2DC1',
            'cf':  'CG2DC1',
            'cg':  'CG1T2',
            'ch':  'CG1T2',
            'cx':  'CG3C31',
            'cy':  'CG3C41',
            'cu':  'CG2D1',
            'cv':  'CG2D1',
            'cz':  'CG2R64',
            'h1':  'HGA1',
            'h2':  'HGA2',
            'h3':  'HGA3',
            'h4':  'HGR62',
            'h5':  'HGA5',
            'ha':  'HGR61',
            'hc':  'HGA3',
            'hn':  'HGP1',
            'ho':  'HGP1',
            'hp':  'HGP5',
            'hs':  'HGP3',
            'hw':  'HW',
            'hx':  'HGPAM1',
            'f':   'FGA1',
            'cl':  'CLGA1',
            'br':  'BRGA1',
            'i':   'IGR1',
            'n':   'NG2S1',
            'n1':  'NG1T1',
            'n2':  'NG2D1',
            'n3':  'NG3N1',
            'n4':  'NG3P1',
            'na':  'NG2R51',
            'nb':  'NG2R60',
            'nc':  'NG2D1',
            'nd':  'NG2R50',
            'nh':  'NG2S3',
            'no':  'NG2O1',
            'ni':  'NG3N1',
            'nj':  'NG3N1',
            'nk':  'NG3P1',
            'nl':  'NG3P1',
            'o':   'OG2D1',
            'oh':  'OG311',
            'os':  'OG301',
            'op':  'OG301',
            'oq':  'OG301',
            'ow':  'OT',
            'p2':  'PG0',
            'p3':  'PG1',
            'p4':  'PG1',
            'p5':  'PG2',
            'pb':  'PG1',
            'pc':  'PG1',
            'pd':  'PG1',
            's':   'SG311',
            's2':  'SG2D1',
            's4':  'SG3O2',
            's6':  'SG3O3',
            'sh':  'SG311',
            'ss':  'SG2R50',
            'sp':  'SG3O1',
            'sq':  'SG3O1'
        }
        
        # Initialize data structures
        self.charmm_bonds: List[CharmmBond] = []
        self.charmm_angles: List[CharmmAngle] = []
        self.charmm_dihedrals: List[CharmmDihedral] = []
        self.charmm_nonbonds: List[CharmmNonbond] = []
        self.frcmod_bonds: List[FrcmodBond] = []
        self.frcmod_angles: List[FrcmodAngle] = []
        self.frcmod_dihedrals: List[FrcmodDihedral] = []
        self.frcmod_nonbonds: List[FrcmodNonbond] = []
        
        # Build lookup dictionaries for faster searching
        self.bond_lookup: Dict[Tuple[str, str], List[CharmmBond]] = {}
        self.angle_lookup: Dict[Tuple[str, str, str], List[CharmmAngle]] = {}
        self.dihedral_lookup: Dict[Tuple[str, str, str, str], List[CharmmDihedral]] = {}
        self.nonbond_lookup: Dict[str, CharmmNonbond] = {}

    #Create a normalized key for bond lookup (order independent)
    def _get_normalized_bond_key(self, atom1: str, atom2: str) -> Tuple[str, str]:
        return tuple(sorted([atom1, atom2]))
    
    #Create a normalized key for angle lookup (direction independent)
    def _get_normalized_angle_key(self, atom1: str, atom2: str, atom3: str) -> Tuple[str, str, str]:
        if atom1 < atom3:
            return (atom1, atom2, atom3)
        return (atom3, atom2, atom1)
    
    #Create a normalized key for dihedral lookup (direction independent)
    def _get_normalized_dihedral_key(self, atom1: str, atom2: str, atom3: str, atom4: str) -> Tuple[str, str, str, str]:
        forward = (atom1, atom2, atom3, atom4)
        reverse = (atom4, atom3, atom2, atom1)
        return forward if forward < reverse else reverse
    
    #Build lookup tables for faster parameter matching
    def _build_lookup_tables(self) -> None:
        # Build bond lookup
        for bond in self.charmm_bonds:
            key = self._get_normalized_bond_key(bond.atom1, bond.atom2)
            if key not in self.bond_lookup:
                self.bond_lookup[key] = []
            self.bond_lookup[key].append(bond)
            
            # Also add prefix-based keys for generic matching
            prefix_key = self._get_normalized_bond_key(bond.atom1[:2], bond.atom2[:2])
            if prefix_key not in self.bond_lookup:
                self.bond_lookup[prefix_key] = []
            self.bond_lookup[prefix_key].append(bond)
        
        # Build angle lookup
        for angle in self.charmm_angles:
            key = self._get_normalized_angle_key(angle.atom1, angle.atom2, angle.atom3)
            if key not in self.angle_lookup:
                self.angle_lookup[key] = []
            self.angle_lookup[key].append(angle)
            
            # Also add prefix-based keys for generic matching
            prefix_key = self._get_normalized_angle_key(
                angle.atom1[:2], angle.atom2[:2], angle.atom3[:2]
            )
            if prefix_key not in self.angle_lookup:
                self.angle_lookup[prefix_key] = []
            self.angle_lookup[prefix_key].append(angle)
        
        # Build dihedral lookup
        for dihedral in self.charmm_dihedrals:
            key = self._get_normalized_dihedral_key(
                dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            )
            if key not in self.dihedral_lookup:
                self.dihedral_lookup[key] = []
            self.dihedral_lookup[key].append(dihedral)
            
            # Also add prefix-based keys for generic matching
            prefix_key = self._get_normalized_dihedral_key(
                dihedral.atom1[:2], dihedral.atom2[:2], dihedral.atom3[:2], dihedral.atom4[:2]
            )
            if prefix_key not in self.dihedral_lookup:
                self.dihedral_lookup[prefix_key] = []
            self.dihedral_lookup[prefix_key].append(dihedral)
        
        # Build nonbond lookup
        for nonbond in self.charmm_nonbonds:
            self.nonbond_lookup[nonbond.atom_type] = nonbond
            # Also add prefix for generic matching
            self.nonbond_lookup[nonbond.atom_type[:2]] = nonbond

    #Read CHARMM parameter file bonds, angles, dihedrals, and nonbonded sections.
    def read_charmm_prm(self, filename: str) -> None:
        
        if not os.path.exists(filename):
            raise FileNotFoundError(f"CHARMM parameter file not found: {filename}")
        
        current_section = None
        line_num = 0
        
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line_num += 1
                    line = line.strip()
                    if not line or line.startswith('!'):
                        continue
                    
                    # Determine section
                    if 'BONDS' in line:
                        current_section = 'BONDS'
                        continue
                    elif 'ANGLES' in line:
                        current_section = 'ANGLES'
                        continue
                    elif 'DIHEDRALS' in line or 'TORSIONS' in line:
                        current_section = 'DIHEDRALS'
                        continue
                    elif 'NONBONDED' in line:
                        current_section = 'NONBONDED'
                        continue
                    elif 'IMPROPER' in line:
                        current_section = None
                        continue
                    
                    # Parse data for current section
                    parts = line.split('!')
                    data = parts[0].split()
                    comment = parts[1].strip() if len(parts) > 1 else ""
                    
                    try:
                        if current_section == 'BONDS' and len(data) >= 4:
                            self.charmm_bonds.append(CharmmBond(
                                atom1=data[0],
                                atom2=data[1],
                                force=float(data[2]),
                                length=float(data[3]),
                                comment=comment
                            ))
                        elif current_section == 'ANGLES' and len(data) >= 5:
                            self.charmm_angles.append(CharmmAngle(
                                atom1=data[0],
                                atom2=data[1],
                                atom3=data[2],
                                force=float(data[3]),
                                angle=float(data[4]),
                                comment=comment
                            ))
                        elif current_section == 'DIHEDRALS' and len(data) >= 7:
                            self.charmm_dihedrals.append(CharmmDihedral(
                                atom1=data[0],
                                atom2=data[1],
                                atom3=data[2],
                                atom4=data[3],
                                force=float(data[4]),
                                periodicity=int(data[5]),
                                phase=float(data[6]),
                                comment=comment
                            ))
                        elif current_section == 'NONBONDED' and len(data) >= 4:
                            if not data[0].startswith('!') and not data[0].startswith('#'):
                                self.charmm_nonbonds.append(CharmmNonbond(
                                    atom_type=data[0],
                                    ignored=float(data[1]),
                                    epsilon=float(data[2]),
                                    rmin_half=float(data[3]),
                                    comment=comment
                                ))
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line {line_num} in {filename}: {line}")
                        logger.warning(f"Exception: {e}")
                        continue
                        
            # Build lookup tables after reading all parameters
            self._build_lookup_tables()
            logger.info(f"Read CHARMM parameter file: {filename}")
            logger.info(f"Found {len(self.charmm_bonds)} bonds, {len(self.charmm_angles)} angles, "
                      f"{len(self.charmm_dihedrals)} dihedrals, {len(self.charmm_nonbonds)} nonbonds")
        
        except Exception as e:
            logger.error(f"Error reading CHARMM parameter file {filename}: {str(e)}")
            raise

    #Read FRCMOD file bonds, angles, dihedrals, and nonbonded sections.
    def read_frcmod_parameters(self, filename: str) -> None:
        
        if not os.path.exists(filename):
            raise FileNotFoundError(f"FRCMOD parameter file not found: {filename}")
        
        current_section = None
        line_num = 0
        
        try:
            with open(filename, 'r') as f:
                for line in f:
                    line_num += 1
                    line = line.strip()
                    if not line or line.startswith('Remark'):
                        continue
                        
                    # Determine section
                    if line.startswith('BOND'):
                        current_section = 'BOND'
                        continue
                    elif line.startswith('ANGLE'):
                        current_section = 'ANGLE'
                        continue
                    elif line.startswith('DIHE'):
                        current_section = 'DIHE'
                        continue
                    elif line.startswith('NONBON'):
                        current_section = 'NONBON'
                        continue
                    elif line.startswith('IMPROPER'):
                        current_section = None  # Skip impropers for now
                        continue
                    
                    try:
                        if current_section == 'BOND':
                            parts = line.split()
                            if len(parts) >= 3:
                                atoms = parts[0].split('-')
                                if len(atoms) == 2:
                                    self.frcmod_bonds.append(FrcmodBond(
                                        atom1=atoms[0],
                                        atom2=atoms[1],
                                        force=float(parts[1]),
                                        length=float(parts[2])
                                    ))
                        
                        elif current_section == 'ANGLE':
                            parts = line.split()
                            if len(parts) >= 3:
                                atoms = parts[0].split('-')
                                if len(atoms) == 3:
                                    self.frcmod_angles.append(FrcmodAngle(
                                        atom1=atoms[0],
                                        atom2=atoms[1],
                                        atom3=atoms[2],
                                        force=float(parts[1]),
                                        angle=float(parts[2])
                                    ))
                        
                        elif current_section == 'DIHE':
                            parts = line.split()
                            if len(parts) >= 4:
                                atoms = parts[0].split('-')
                                if len(atoms) == 4:
                                    self.frcmod_dihedrals.append(FrcmodDihedral(
                                        atom1=atoms[0],
                                        atom2=atoms[1],
                                        atom3=atoms[2],
                                        atom4=atoms[3],
                                        periodicity=int(parts[1]),
                                        force=float(parts[2]),
                                        phase=float(parts[3]),
                                        divider=float(parts[4]) if len(parts) > 4 else 1.0
                                    ))
                        
                        elif current_section == 'NONBON':
                            parts = line.split()
                            if len(parts) >= 3:
                                self.frcmod_nonbonds.append(FrcmodNonbond(
                                    atom_type=parts[0],
                                    radius=float(parts[1]),
                                    well_depth=float(parts[2])
                                ))
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line {line_num} in {filename}: {line}")
                        logger.warning(f"Exception: {e}")
                        continue
            
            logger.info(f"Read FRCMOD parameter file: {filename}")
            logger.info(f"Found {len(self.frcmod_bonds)} bonds, {len(self.frcmod_angles)} angles, "
                      f"{len(self.frcmod_dihedrals)} dihedrals, {len(self.frcmod_nonbonds)} nonbonds")
            
            # Track missing atom types
            self._check_amber_atom_types()
        
        except Exception as e:
            logger.error(f"Error reading FRCMOD parameter file {filename}: {str(e)}")
            raise

    #Check for AMBER atom types not in the conversion dictionary
    def _check_amber_atom_types(self) -> None:
        found_types = set()
        
        # Collect all unique AMBER atom types from the frcmod parameters
        for bond in self.frcmod_bonds:
            found_types.add(bond.atom1)
            found_types.add(bond.atom2)
        
        for angle in self.frcmod_angles:
            found_types.add(angle.atom1)
            found_types.add(angle.atom2)
            found_types.add(angle.atom3)
            
        for dihedral in self.frcmod_dihedrals:
            found_types.add(dihedral.atom1)
            found_types.add(dihedral.atom2)
            found_types.add(dihedral.atom3)
            found_types.add(dihedral.atom4)
            
        for nonbond in self.frcmod_nonbonds:
            found_types.add(nonbond.atom_type)
        
        # Check which types are missing from our conversion dictionary
        missing_types = [t for t in found_types if t not in self.amber_to_charmm]
        
        if missing_types:
            logger.warning(f"Found {len(missing_types)} AMBER atom types without CHARMM equivalents:")
            logger.warning(", ".join(missing_types))

    #Find matching bond in CHARMM parameters.
    def find_matching_bond(self, type1: str, type2: str, target_length: float) -> Optional[CharmmBond]:
        
        charmm_type1 = self.amber_to_charmm.get(type1)
        charmm_type2 = self.amber_to_charmm.get(type2)
        
        if not (charmm_type1 and charmm_type2):
            return None
        
        # Try exact match
        key = self._get_normalized_bond_key(charmm_type1, charmm_type2)
        if key in self.bond_lookup:
            # Find the best match based on bond length
            best_match = None
            min_diff = float('inf')
            
            for bond in self.bond_lookup[key]:
                diff = abs(bond.length - target_length)
                if diff < min_diff:
                    min_diff = diff
                    best_match = bond
                    
            if best_match and min_diff < 0.1:
                return best_match
        
        # Try generic match (using prefixes)
        prefix_key = self._get_normalized_bond_key(charmm_type1[:2], charmm_type2[:2])
        if prefix_key in self.bond_lookup:
            best_match = None
            min_diff = float('inf')
            
            for bond in self.bond_lookup[prefix_key]:
                diff = abs(bond.length - target_length)
                if diff < min_diff:
                    min_diff = diff
                    best_match = bond
                    
            if best_match and min_diff < 0.1:
                return best_match
                
        return None
    
    #Find matching angle in CHARMM parameters.
    def find_matching_angle(self, type1: str, type2: str, type3: str, target_angle: float) -> Optional[CharmmAngle]:
        
        charmm_type1 = self.amber_to_charmm.get(type1)
        charmm_type2 = self.amber_to_charmm.get(type2)
        charmm_type3 = self.amber_to_charmm.get(type3)
    
        # If any atom type doesn't have a CHARMM equivalent, return None
        if not (charmm_type1 and charmm_type2 and charmm_type3):
            return None
        
        # Try exact match
        key = self._get_normalized_angle_key(charmm_type1, charmm_type2, charmm_type3)
        if key in self.angle_lookup:
            # Find the best match based on angle value
            best_match = None
            min_diff = float('inf')
            
            for angle in self.angle_lookup[key]:
                diff = abs(angle.angle - target_angle)
                if diff < min_diff:
                    min_diff = diff
                    best_match = angle
                    
            if best_match and min_diff <= 2.0:
                return best_match
        
        # Try generic match (using prefixes)
        prefix_key = self._get_normalized_angle_key(
            charmm_type1[:2], charmm_type2[:2], charmm_type3[:2]
        )
        if prefix_key in self.angle_lookup:
            best_match = None
            min_diff = float('inf')
            
            for angle in self.angle_lookup[prefix_key]:
                diff = abs(angle.angle - target_angle)
                if diff < min_diff:
                    min_diff = diff
                    best_match = angle
                    
            if best_match and min_diff <= 2.0:
                return best_match
            
        return None

    #Find matching dihedral in CHARMM parameters.
    def find_matching_dihedral(self, type1: str, type2: str, type3: str, type4: str) -> Optional[CharmmDihedral]:
        
        charmm_type1 = self.amber_to_charmm.get(type1)
        charmm_type2 = self.amber_to_charmm.get(type2)
        charmm_type3 = self.amber_to_charmm.get(type3)
        charmm_type4 = self.amber_to_charmm.get(type4)
        
        # If any atom type doesn't have a CHARMM equivalent, return None
        if not all([charmm_type1, charmm_type2, charmm_type3, charmm_type4]):
            return None
        
        # Try exact match
        key = self._get_normalized_dihedral_key(
            charmm_type1, charmm_type2, charmm_type3, charmm_type4
        )
        if key in self.dihedral_lookup:
            return self.dihedral_lookup[key][0]  # Return first match
        
        # Try generic match (using prefixes)
        prefix_key = self._get_normalized_dihedral_key(
            charmm_type1[:2], charmm_type2[:2], charmm_type3[:2], charmm_type4[:2]
        )
        if prefix_key in self.dihedral_lookup:
            return self.dihedral_lookup[prefix_key][0]  # Return first match
        
        return None

    #Find matching nonbond parameters in CHARMM parameters.
    def find_matching_nonbond(self, atom_type: str) -> Optional[CharmmNonbond]:
        
        # First, determine the proper CHARMM atom type
        if atom_type in self.amber_to_charmm:
            charmm_type = self.amber_to_charmm[atom_type]
        elif atom_type[0].islower() and any(c.isdigit() for c in atom_type):
            # Strip numbers and look up the base type
            base_type = ''.join(c for c in atom_type if not c.isdigit()).lower()
            charmm_type = self.amber_to_charmm.get(base_type)
        else:
            # No conversion for atom types like metal that don't have a mapping
            charmm_type = None
        
        # Search for the exact CHARMM type
        if charmm_type and charmm_type in self.nonbond_lookup:
            return self.nonbond_lookup[charmm_type]
        
        # If no exact match, try a prefix match
        if charmm_type and charmm_type[:2] in self.nonbond_lookup:
            return self.nonbond_lookup[charmm_type[:2]]
        
        return None

    #Convert a nonbond parameter to CHARMM format.
    def convert_nonbond_parameters(self, nonbond: FrcmodNonbond) -> str:
        
        charmm_nonbond = self.find_matching_nonbond(nonbond.atom_type)
        
        if charmm_nonbond:
            # Found a matching CHARMM parameter
            return f"{nonbond.atom_type:<10} {charmm_nonbond.ignored:>8.3f} {charmm_nonbond.epsilon:>10.4f} {charmm_nonbond.rmin_half:>10.4f}"
        else:
            # No match found - handle according to rules
            if any(c.isdigit() for c in nonbond.atom_type):
                # Has a number - remove it and try again
                base_type = ''.join(c for c in nonbond.atom_type if not c.isdigit()).lower()
                charmm_nonbond = self.find_matching_nonbond(base_type)
                if charmm_nonbond:
                    return f"{nonbond.atom_type:<10} {charmm_nonbond.ignored:>8.3f} {charmm_nonbond.epsilon:>10.4f} {charmm_nonbond.rmin_half:>10.4f}"
            
            # Default case - transform directly from frcmod format
            # Negate well_depth for CHARMM format
            return f"{nonbond.atom_type:<10} {0.0:>8.3f} {-nonbond.well_depth:>10.4f} {nonbond.radius:>10.4f}"

    #Convert AMBER parameters to CHARMM format.
    def convert_parameters(self) -> str:
        
        output = []
        converted_count = {'bonds': 0, 'angles': 0, 'dihedrals': 0, 'nonbonds': 0}
        missed_count = {'bonds': 0, 'angles': 0, 'dihedrals': 0, 'nonbonds': 0}
    
        # Convert bonds
        output.append("BONDS")
        for bond in self.frcmod_bonds:
            charmm_bond = self.find_matching_bond(bond.atom1, bond.atom2, bond.length)

            if charmm_bond:
                output.append(f"{bond.atom1:<4} {bond.atom2:<4} {charmm_bond.force:>8.2f} {bond.length:>8.3f}")
            else:
                output.append(f"{bond.atom1:<4} {bond.atom2:<4} {bond.force:>8.2f} {bond.length:>8.3f}")
    
        # Convert angles
        output.append("\nANGLES")
        for angle in self.frcmod_angles:
            charmm_angle = self.find_matching_angle(angle.atom1, angle.atom2, angle.atom3, angle.angle)
        
            if charmm_angle:
                # Use CHARMM force constant with original equilibrium angle
                output.append(f"{angle.atom1:<4} {angle.atom2:<4} {angle.atom3:<4} {charmm_angle.force:>8.2f} {angle.angle:>8.3f}")
            else:
                # If no CHARMM types available or no match within 2.0 degrees, use half of original force
                output.append(f"{angle.atom1:<4} {angle.atom2:<4} {angle.atom3:<4} {angle.force/2:>8.2f} {angle.angle:>8.3f}")
        
        # Convert dihedrals
        output.append("\nDIHEDRALS")
        for dihedral in self.frcmod_dihedrals:
            charmm_dihedral = self.find_matching_dihedral(
                dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4
            )
            
            if charmm_dihedral and all(atom in self.amber_to_charmm for atom in 
                                      [dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4]):
                # Use CHARMM force constant with original phase
                output.append(f"{dihedral.atom1:<4} {dihedral.atom2:<4} {dihedral.atom3:<4} {dihedral.atom4:<4} "
                             f"{charmm_dihedral.force:>8.4f} {dihedral.periodicity:>3d} {dihedral.phase:>8.4f}")
            else:
                # If no CHARMM match, divide force by 2 and use periodicity of 2
                force = dihedral.force / 2.0
                output.append(f"{dihedral.atom1:<4} {dihedral.atom2:<4} {dihedral.atom3:<4} {dihedral.atom4:<4} "
                             f"{force:>8.4f} {dihedral.periodicity:>3d} {dihedral.phase:>8.4f}")
    
        # Convert nonbonds
        output.append("\nNONBONDED")
        for nonbond in self.frcmod_nonbonds:
            output.append(self.convert_nonbond_parameters(nonbond))
    
        # Add a blank line and END marker
        output.append("\n")
    
        return "\n".join(output)

def main():
    converter = ParameterConverter()
   
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    charmm_file_path = os.path.join(script_dir, "charmm.prm")

    # Read input files
    converter.read_charmm_prm(charmm_file_path)

    converter.read_frcmod_parameters("COMPLEX.frcmod")
    
    # Convert and write output
    output = converter.convert_parameters()
    with open("COMPLEX.prm", "w") as f:
        f.write(output)

if __name__ == "__main__":
    main()
