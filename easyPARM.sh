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
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


# Colors
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Print without string interpolation
echo -e "${PURPLE}════════════════════════════════════════════════════════════════════════${NC}"
echo -e "${CYAN}
███████╗ █████╗ ███████╗██╗   ██╗██████╗  █████╗ ██████╗ ███╗   ███╗
██╔════╝██╔══██╗██╔════╝╚██╗ ██╔╝██╔══██╗██╔══██╗██╔══██╗████╗ ████║
███████╗███████║███████╗ ╚████╔╝ ██████╔╝███████║██████╔╝██╔████╔██║
██╔══╝  ██╔══██║╚════██║  ╚██╔╝  ██╔═══╝ ██╔══██║██╔══██╗██║╚██╔╝██║
███████╗██║  ██║███████║   ██║   ██║     ██║  ██║██║  ██║██║ ╚═╝ ██║
╚══════╝╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝${NC}"
echo -e "${GREEN}                ⚡ Automated Force Fields for Metals ⚡${NC}"
echo -e "${PURPLE}════════════════════════════════════════════════════════════════════════${NC}"
echo -e "                     ${YELLOW}Version 4.25 — May 2026${NC}"
echo -e "${PURPLE}════════════════════════════════════════════════════════════════════════${NC}"


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

# Capture the directory from which this script is called
ORIGINAL_DIR="$(pwd)"

# Determine the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate to the scripts directory
cd "$SCRIPT_DIR/scripts"

# Check if an input file was provided with .inp or .input extension
if [ $# -eq 1 ] && [[ "$1" == *.inp || "$1" == *.input || "$1" == *.config ]]; then
    # Non-interactive mode
    input_file="$1"

    # Use Python script to convert the input file to the format expected by 01_easyPARM.sh
    temp_input="tempx.input"
    rm -f "$ORIGINAL_DIR/$temp_input"
    echo " "
    python3 "$SCRIPT_DIR/scripts/input_generator.py" "$ORIGINAL_DIR/$input_file" "$ORIGINAL_DIR/$temp_input"

    # Run 01_easyPARM.sh with the converted input
    if [ -f "$ORIGINAL_DIR/$temp_input" ]; then
	    echo -e "\n\033[1;32m=== NON-INTERACTIVE MODE ACTIVATED ===\033[0m"
	    echo -e "\033[1;90mRunning easyPARM with configuration from $input_file\033[0m\n"
	    ./01_easyPARM.sh "$ORIGINAL_DIR" < "$ORIGINAL_DIR/$temp_input"
	    rm -f "$ORIGINAL_DIR/$temp_input"
    else
	    # Interactive mode 
	    echo -e "\n\033[1;34m=== INTERACTIVE MODE ACTIVATED ===\033[0m"
	    echo -e "\033[1;90mRunning easyPARM with user prompts\033[0m\n"
	    ./01_easyPARM.sh "$ORIGINAL_DIR"
    fi 
else
    # Interactive mod
    echo -e "\n\033[1;34m=== INTERACTIVE MODE ACTIVATED ===\033[0m"
    echo -e "\033[1;90mRunning easyPARM with user prompts\033[0m\n"
    echo " "
    echo "================================="
    echo "          easyPARM Menu          "
    echo "================================="
    echo "Select your option:"
    echo "1- Generate molecular complex parameters"
    echo "2- Generate metalloprotein .xyz structure"
    echo "3- Generate metal–nucleic acid system .xyz structure"
    echo "4- Convert AMBER parameters to OpenMM or GROMACS format"
    choice=$(get_valid_input "Enter your choice: " "1 2 3 4")
    if [ "$choice" = "1" ]; then
        ./01_easyPARM.sh "$ORIGINAL_DIR"
    elif [ "$choice" = "2" ]; then
        ./extract_metal_coordination.sh "$ORIGINAL_DIR" 1
    elif [ "$choice" = "3" ]; then
        ./extract_metal_coordination.sh "$ORIGINAL_DIR" 2
    elif [ "$choice" = "4" ]; then
        ./amber_converter.sh "$ORIGINAL_DIR"
    fi
fi
