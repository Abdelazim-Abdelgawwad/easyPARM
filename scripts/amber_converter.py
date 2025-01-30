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


import parmed as pmd
import openmm
from openmm import app
from openmm import unit
import sys
import os

#Convert AMBER files to OpenMM system with optional periodic boundary conditions.
def convert_amber_to_openmm(prmtop_file, inpcrd_file, add_box_if_absent=True):
    # Load AMBER files using ParmEd
    amber_struct = pmd.load_file(prmtop_file, inpcrd_file)

    # Check if structure has a periodic box
    has_box = amber_struct.box is not None

    if not has_box and add_box_if_absent:
        # Add a default cubic periodic box based on the system size
        max_dimension = max(amber_struct.coordinates.max(axis=0) - amber_struct.coordinates.min(axis=0))
        amber_struct.box = [max_dimension + 10, max_dimension + 10, max_dimension + 10, 90.0, 90.0, 90.0]  # Box dimensions in Å

    # Choose appropriate nonbonded method based on periodicity
    if amber_struct.box is not None:
        nonbonded_method = app.PME
    else:
        nonbonded_method = app.NoCutoff

    # Convert to OpenMM System
    system = amber_struct.createSystem(
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=1.0 * unit.nanometers if amber_struct.box is not None else None,
        constraints=app.HBonds
    )

    # Get OpenMM Topology
    topology = amber_struct.topology

    # Get positions (coordinates)
    positions = amber_struct.positions

    return system, topology, positions

#Save OpenMM system to XML file.
def save_openmm_xml(system, output_file):
    with open(output_file, 'w') as f:
        xml = openmm.XmlSerializer.serialize(system)
        f.write(xml)

#Convert AMBER parameter and coordinate files to GROMACS format.
def convert_amber_to_gromacs(amber_prmtop, amber_inpcrd, output_prefix):
    try:
        # Load the AMBER files
        amber = pmd.load_file(amber_prmtop, amber_inpcrd)

        # Create output filenames
        gro_file = f"{output_prefix}.gro"
        top_file = f"{output_prefix}.top"

        # Save as GROMACS format
        amber.save(gro_file, format='gro')
        amber.save(top_file, format='gromacs')

    except Exception as e:
        print(f"Error during conversion: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(1)

    prmtop_file = sys.argv[1]
    inpcrd_file = sys.argv[2]
    mode = int(sys.argv[3])

    if mode == 1:
        output_xml = "system.xml"

        # Convert files, preserving the original box size if present
        system, topology, positions = convert_amber_to_openmm(
            prmtop_file,
            inpcrd_file,
            add_box_if_absent=True  # Add a periodic box only if it is missing
        )

        # Save system to XML
        save_openmm_xml(system, output_xml)

    elif mode == 2:
        output_prefix = "system_gmx"

        # Convert files to GROMACS format
        convert_amber_to_gromacs(prmtop_file, inpcrd_file, output_prefix)

    else:
        print("Invalid mode. Use 1 for AMBER to OpenMM or 2 for AMBER to GROMACS.")

