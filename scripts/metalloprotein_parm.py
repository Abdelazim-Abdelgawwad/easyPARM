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
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import re
import os
import sys


# Read metalloprotein_atomtype.dat file and return mapping from new_atom_type to original_atom_type
def read_atom_type_mapping(filename):
    new_to_original = {}
    original_to_new = {}

    try:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    new_atom_type = parts[0]
                    original_atom_type = parts[1]
                    new_to_original[new_atom_type] = original_atom_type
                    original_to_new[original_atom_type] = new_atom_type
    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
        exit(1)

    return new_to_original, original_to_new


# Read protein_parm.dat file and extract bond and angle parameters
def read_protein_parm_data(filename):
    bond_params = {}
    angle_params = {}

    try:
        with open(filename, 'r') as f:
            section = None

            for line in f:
                line = line.strip()

                if line.startswith('BOND'):
                    section = 'BOND'
                    continue
                elif line.startswith('ANGLE'):
                    section = 'ANGLE'
                    continue
                elif line.startswith('DIHE'):
                    section = 'DIHE'
                    continue

                if not line or line.startswith('#'):
                    continue

                if section == 'BOND':
                    # \s*-\s* handles both "C-O" and "C -O" formats
                    match = re.match(r'(\S+)\s*-\s*(\S+)\s+(\S+)\s+(\S+)', line)
                    if match:
                        atom1, atom2, force, dist = match.groups()
                        key1 = f"{atom1}-{atom2}"
                        key2 = f"{atom2}-{atom1}"
                        bond_params[key1] = (float(force), float(dist))
                        bond_params[key2] = (float(force), float(dist))

                elif section == 'ANGLE':
                    # \s*-\s* handles both "C-O-CT" and "C -O -CT" formats
                    match = re.match(r'(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s+(\S+)\s+(\S+)', line)
                    if match:
                        atom1, atom2, atom3, force, angle = match.groups()
                        key1 = f"{atom1}-{atom2}-{atom3}"
                        key2 = f"{atom3}-{atom2}-{atom1}"
                        angle_params[key1] = (float(force), float(angle))
                        angle_params[key2] = (float(force), float(angle))

    except FileNotFoundError:
        print(f"Error: Could not find file {filename}")
        exit(1)

    return bond_params, angle_params


def replace_numbers_in_line(original_line, new_values):
    """
    Replace only the numeric parameter values in a line, preserving
    all atom types, spacing, and formatting exactly as they are.

    Finds all numeric tokens (integers or floats) after the atom-type
    string and replaces them one-by-one with the new values.

    Example:
        original_line = "C -O2   622.90  1.225\n"
        new_values    = [570.0, 1.229]
        result        = "C -O2   570.00  1.229\n"
    """
    # Split line into: everything up to first digit block, then numeric tokens
    # We locate each number by its span and replace in reverse order to keep positions valid
    number_pattern = re.compile(r'\d+\.\d+|\d+')

    # Find the atom-type prefix: everything up to the first whitespace-separated number
    # Strategy: find all number matches, replace them with new_values in order
    matches = list(number_pattern.finditer(original_line))

    if len(matches) < len(new_values):
        # Not enough numbers found  return original unchanged
        return original_line

    # Build the new line by replacing numbers from right to left
    # (so earlier positions are not shifted)
    result = original_line
    # Only replace as many numbers as we have new values, taken from the end
    # of the type string  i.e. the last len(new_values) matches
    targets = matches[-len(new_values):]

    for match, new_val in zip(reversed(targets), reversed(new_values)):
        original_num = match.group()
        # Format the new value with the same number of decimal places as original
        if '.' in original_num:
            decimal_places = len(original_num.split('.')[1])
            formatted = f"{new_val:.{decimal_places}f}"
        else:
            formatted = str(int(new_val))

        # Pad/trim to same width as original to preserve column alignment
        if len(formatted) < len(original_num):
            formatted = formatted.rjust(len(original_num))

        result = result[:match.start()] + formatted + result[match.end():]

    return result


# Update COMPLEX.frcmod file with parameters from protein_parm.dat based on atom type mappings.
def update_complex_frcmod(frcmod_filename, new_to_original, original_to_new, protein_bond_params, protein_angle_params):
    try:
        with open(frcmod_filename, 'r') as f:
            frcmod_content = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find file {frcmod_filename}")
        exit(1)

    section = None
    updated_content = []

    for line in frcmod_content:
        original_line = line
        stripped_line = line.strip()

        # Section detection
        if stripped_line == 'BOND':
            section = 'BOND'
            updated_content.append(original_line)
            continue
        elif stripped_line == 'ANGLE':
            section = 'ANGLE'
            updated_content.append(original_line)
            continue
        elif stripped_line == 'DIHE':
            section = 'DIHE'
            updated_content.append(original_line)
            continue
        elif not stripped_line:
            updated_content.append(original_line)
            continue

        if section == 'BOND':
            match = re.match(r'(\S+)\s*-\s*(\S+)\s+(\S+)\s+(\S+)', stripped_line)
            if match:
                atom1, atom2, force, dist = match.groups()

                # Map new atom types ? original atom types for lookup
                orig_atom1 = new_to_original.get(atom1, atom1)
                orig_atom2 = new_to_original.get(atom2, atom2)

                protein_key1 = f"{orig_atom1}-{orig_atom2}"
                protein_key2 = f"{orig_atom2}-{orig_atom1}"

                if protein_key1 in protein_bond_params or protein_key2 in protein_bond_params:
                    key = protein_key1 if protein_key1 in protein_bond_params else protein_key2
                    new_force, new_dist = protein_bond_params[key]
                    # Replace only the numbers  keep atom types and spacing identical
                    new_line = replace_numbers_in_line(original_line, [new_force, new_dist])
                    updated_content.append(new_line)
                else:
                    updated_content.append(original_line)
            else:
                print(f"WARNING: Could not parse BOND line: {repr(stripped_line)}")
                updated_content.append(original_line)

        elif section == 'ANGLE':
            match = re.match(r'(\S+)\s*-\s*(\S+)\s*-\s*(\S+)\s+(\S+)\s+(\S+)', stripped_line)
            if match:
                atom1, atom2, atom3, force, angle = match.groups()

                # Map new atom types ? original atom types for lookup
                orig_atom1 = new_to_original.get(atom1, atom1)
                orig_atom2 = new_to_original.get(atom2, atom2)
                orig_atom3 = new_to_original.get(atom3, atom3)

                protein_key1 = f"{orig_atom1}-{orig_atom2}-{orig_atom3}"
                protein_key2 = f"{orig_atom3}-{orig_atom2}-{orig_atom1}"

                if protein_key1 in protein_angle_params or protein_key2 in protein_angle_params:
                    key = protein_key1 if protein_key1 in protein_angle_params else protein_key2
                    new_force, new_angle = protein_angle_params[key]
                    # Replace only the numbers  keep atom types and spacing identical
                    new_line = replace_numbers_in_line(original_line, [new_force, new_angle])
                    updated_content.append(new_line)
                else:
                    updated_content.append(original_line)
            else:
                print(f"WARNING: Could not parse ANGLE line: {repr(stripped_line)}")
                updated_content.append(original_line)

        else:
            updated_content.append(original_line)

    # Write backup then updated file
    backup_filename = f"{frcmod_filename}.bak"
    try:
        with open(backup_filename, 'w') as f:
            f.writelines(frcmod_content)

        with open(frcmod_filename, 'w') as f:
            f.writelines(updated_content)

    except Exception as e:
        print(f"Error updating file: {e}")
        exit(1)


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 code.py <parm_file>")
        exit(1)

    parm_type = sys.argv[1]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parm_file = os.path.join(script_dir, parm_type)

    if not os.path.exists(parm_type):
        print(f"Error: Parm file '{parm_type}' not found")
        sys.exit(1)

    atom_type_file = "metalloprotein_atomtype.dat"
    frcmod_file = "COMPLEX.frcmod"

    new_to_original, original_to_new = read_atom_type_mapping(atom_type_file)
    protein_bond_params, protein_angle_params = read_protein_parm_data(parm_file)
    update_complex_frcmod(frcmod_file, new_to_original, original_to_new, protein_bond_params, protein_angle_params)


if __name__ == "__main__":
    main()

