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
#                              |  $$$$$$/              Ver. 3.10 - 12 February 2025                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

echo "  "

# Detect the directory where the current script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the directory from which easyPARM.sh was executed 
RUN_DIR="$1"


cd "$RUN_DIR"


echo "Running in directory: $RUN_DIR"
echo "Script directory: $SCRIPT_DIR"
echo " "

read -p "Please provide the metalloprotein pdb file: " protein_pdb

if [ -f "$SCRIPT_DIR/prepare_metalloprotein_xyz.py" ]; then
        python3 "$SCRIPT_DIR/prepare_metalloprotein_xyz.py" "$RUN_DIR/$protein_pdb"
    else
        echo "Script prepare_metalloprotein_xyz.py not found in $SCRIPT_DIR. Exiting."
        exit 1
fi

    # Check if the script was successful
if [ $? -ne 0 ]; then
        echo "Failed to execute prepare_metalloprotein_xyz.py. Exiting."
        exit 1
fi

echo ""
echo "XYZ Output: initial_structure.xyz"
echo ""
echo "You can now proceed to the next step, such as optimization and frequency calculations, to generate the required files for easyPARM."
