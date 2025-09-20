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
#                              |  $$$$$$/              Ver. 4.10 - 20 September 2025                              #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


import sys
import numpy as np
import periodictable

#Extract atomic coordinates, atomic numbers, and the number of atoms
def parse_coordinates(file_path):
    coordinates = []
    atom_names = []
    reading_atoms = False
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() == '$atoms':
                reading_atoms = True
                num_atoms = int(next(file).strip())
                continue
                
            if reading_atoms:
                if line.strip().startswith('$'):
                    break
                if not line.strip():
                    continue
                    
                parts = line.split()
                if len(parts) >= 5:
                    atom_names.append(parts[0].strip())
                    coordinates.append([float(x) for x in parts[-3:]])
    
    coordinates = np.array(coordinates)
    atomic_numbers = np.array([periodictable.elements.symbol(name).number for name in atom_names])
    
    return coordinates, atomic_numbers, atom_names

#Read similar.dat and return a mapping of row_number to zeros_value
def parse_similar_dat(file_path):
    mapping = {}
    with open(file_path, 'r') as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2:
                try:
                    row_number = int(tokens[0])
                    zeros_value = int(tokens[1])
                    mapping[row_number] = zeros_value
                except ValueError:
                    continue
    return mapping

#Write the RESP input file with the given charge and atomic numbers.
def write_resp_input(file_path, charge, atomic_numbers, zeros_mapping):
    with open(file_path, 'w') as f:
        f.write("TITLE\n")
        f.write(" &cntrl\n")
        f.write(" nmol=1,\n")
        f.write(" ihfree=1,\n")
        f.write(" ioutopt = 1,\n")
        f.write(" qwt=0.00050,\n")
        f.write(" /\n")
        f.write("    1.0\n")
        f.write("TITLE\n")
        f.write(f"   {charge}   {len(atomic_numbers)}\n")
        f.write(f"   {charge}   {len(atomic_numbers)}\n")

        for i, atomic_number in enumerate(atomic_numbers, start=1):
            zeros_value = zeros_mapping.get(i, 0)
            f.write(f"   {atomic_number}    {zeros_value}\n")
        f.write("\n\n\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 resp.py <file_path_hessian> <charge> <similar_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    charge = int(sys.argv[2])
    similar_file = sys.argv[3]

    resp_file = "resp.in"

    coordinates, atomic_numbers, atom_names = parse_coordinates(file_path)
    zeros_mapping = parse_similar_dat(similar_file)

    write_resp_input(resp_file, charge, atomic_numbers, zeros_mapping)


