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
#                              |  $$$$$$/              Ver. 3.20 - 3 April 2025                                   #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


import sys

#Extract atomic coordinates, atomic numbers, and the number of atoms
def parse_coordinates(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    start_index = None
    end_index = None
    for i, line in enumerate(lines):
        if '$DATA' in line:
            start_index = i
        elif '$END' in line and start_index is not None:
            end_index = i
            break

    if start_index is None or end_index is None:
        raise ValueError("Coordinates section not found in the file.")

    data_section = lines[start_index + 2:end_index]

    coordinates = []
    atomic_numbers = []

    for line in data_section:
        tokens = line.split()
        if len(tokens) >= 5:
            try:
                atomic_number = int(float(tokens[1]))
                x, y, z = map(float, tokens[2:5])
                coordinates.append((x, y, z))
                atomic_numbers.append(atomic_number)
            except (ValueError, IndexError):
                continue

    num_atoms = len(atomic_numbers)
    if num_atoms == 0:
        raise ValueError("No atomic coordinates found in the $DATA section.")

    return atomic_numbers, coordinates, num_atoms

#Extract the electrostatic potential data.
def parse_electrostatic_potential(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    start_index = None
    for i, line in enumerate(lines):
        if "ELECTROSTATIC POTENTIAL," in line:
            start_index = i + 2
            break

    if start_index is None:
        raise ValueError("Electrostatic potential section not found in the file.")

    potential_section = []
    for line in lines[start_index:]:
        if line.strip() == "" or "$END" in line:
            break
        potential_section.append(line)

    potential_data = []
    for line in potential_section:
        tokens = line.split()
        if len(tokens) >= 5:
            try:
                ipt = int(tokens[0])
                x, y, z = map(float, tokens[1:4])
                elpott = float(tokens[4])
                potential_data.append((ipt, x, y, z, elpott))
            except (ValueError, IndexError):
                continue

    if not potential_data:
        raise ValueError("No electrostatic potential data found.")

    return potential_data

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

#Convert coordinates from Angstrom to Bohr
def convert_to_bohr(coordinates):
    angstrom_to_bohr = 1 / 0.529177249
    return [(x * angstrom_to_bohr, y * angstrom_to_bohr, z * angstrom_to_bohr) for x, y, z in coordinates]

#Write the atomic coordinates and electrostatic potential data to a file
def write_output(file_path, atomic_numbers, coordinates_bohr, potential_data):
    with open(file_path, 'w') as f:
        f.write(f"   {len(atomic_numbers)}{len(potential_data):4d}\n")

        for x, y, z in coordinates_bohr:
            f.write(f"                   {x: .6E}   {y: .6E}   {z: .6E}\n")

        for _, x, y, z, elpott in potential_data:
            x_bohr, y_bohr, z_bohr = x, y, z 
            f.write(f"   {elpott: .6E}   {x_bohr: .6E}   {y_bohr: .6E}   {z_bohr: .6E}\n")

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
        f.write(f"    {charge}   {len(atomic_numbers)}\n")

        for i, atomic_number in enumerate(atomic_numbers, start=1):
            zeros_value = zeros_mapping.get(i, 0)
            f.write(f"   {atomic_number}    {zeros_value}\n")
        f.write("\n\n\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 resp.py <file_path> <charge> <similar_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    charge = int(sys.argv[2])
    similar_file = sys.argv[3]

    output_file = "esp.in"
    resp_file = "resp.in"

    atomic_numbers, coordinates, num_atoms = parse_coordinates(file_path)
    coordinates_bohr = convert_to_bohr(coordinates)
    potential_data = parse_electrostatic_potential(file_path)
    zeros_mapping = parse_similar_dat(similar_file)

    write_output(output_file, atomic_numbers, coordinates_bohr, potential_data)
    write_resp_input(resp_file, charge, atomic_numbers, zeros_mapping)


