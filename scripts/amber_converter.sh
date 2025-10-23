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


get_valid_input() {
    local prompt="$1"
    local valid_options="$2"
    local user_input
    while true; do
        read -p "$prompt" user_input
        # Trim any leading or trailing whitespace
        user_input=$(echo "$user_input" | xargs)
        
        # Check if user_input is within valid_options (by checking exact match)
        if [[ " $valid_options " =~ " $user_input " ]]; then
            echo "$user_input"
            return
        else
            echo "Invalid input. Please try again." >&2
        fi
    done
}


# Detect the directory where the current script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the directory from which easyPARM.sh was executed 
RUN_DIR="$1"


cd "$RUN_DIR"


echo "  "
echo "Running in directory: $RUN_DIR"
echo "Script directory: $SCRIPT_DIR"
echo " "

echo "================================="
echo "      AMBER Converter Menu       "
echo "================================="
echo "Select your option:"
echo "1- AMBER to OpenMM"
echo "2- AMBER to GROMACS"

choice=$(get_valid_input "Enter your choice: " "1 2")
echo " "
read -p "Please provide the prmtop file: " prmtop_file
echo " "
read -p "Please provide the inpcrd file: " inpcrd_file

if [ "$choice" = "1" ]; then
	if [ -f "$SCRIPT_DIR/amber_converter.py" ]; then
        	python3 "$SCRIPT_DIR/amber_converter.py" "$RUN_DIR/$prmtop_file" "$RUN_DIR/$inpcrd_file" 1 > "$RUN_DIR/temp.dat"
		rm "$RUN_DIR/temp.dat"
		echo " "
		echo "The OPENMM parameter is as follows: system.xml"
    	else
        	echo "Script amber_converter.py not found in $SCRIPT_DIR. Exiting."
        	exit 1
	fi

    # Check if the script was successful
	if [ $? -ne 0 ]; then
        	echo "Failed to execute amber_converter.py Exiting."
        	exit 1
	fi

elif [ "$choice" = "2" ]; then
	if [ -f "$SCRIPT_DIR/amber_converter.py" ]; then
        	python3 "$SCRIPT_DIR/amber_converter.py" "$RUN_DIR/$prmtop_file" "$RUN_DIR/$inpcrd_file" 2 > "$RUN_DIR/temp.dat"
		rm "$RUN_DIR/temp.dat"
		echo " "
		echo "The GROMACS parameters are as follows: 1- system_gmx.gro 2- system_gmx.top"
    	else
        	echo "Script amber_converter.py not found in $SCRIPT_DIR. Exiting."
        	exit 1
	fi

    # Check if the script was successful
	if [ $? -ne 0 ]; then
        	echo "Failed to execute amber_converter.py. Exiting."
        	exit 1
	fi

fi   

