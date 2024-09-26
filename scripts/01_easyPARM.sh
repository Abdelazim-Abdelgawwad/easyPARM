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
#                              |  $$$$$$/              Ver. 1.00 - 26 September 2024                              #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, ValÃ¨ncia 46071, Spain          #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de ValÃ¨ncia. E-mail: abdelazim.abdelgawwad@uv.es        #
###################################################################################################################

# Detect the directory where the current script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the directory from which easyPARM.sh was executed 
RUN_DIR="$1"

echo "Running in directory: $RUN_DIR"
echo "Script directory: $SCRIPT_DIR"

# Ensure all outputs go to RUN_DIR
cd "$RUN_DIR"

# Function to ask the user for input again if the information provided was incorrect

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

# Handle user's choice for Amber
echo "Select your option:"
echo "1- Amber is already loaded"
echo "2- Provide the Amber path"
echo " "

while true; do
    choice=$(get_valid_input "Enter your choice (1 or 2): " "1 2")
    
    if [ "$choice" = "1" ]; then
        echo " "
        echo "Amber is assumed to be already loaded. Skipping sourcing."
        break
    elif [ "$choice" = "2" ]; then
        # Prompt the user for the path to Amber
        echo " "
        read -p "Please provide the path for Amber: " amber_path
        
        # Check if the provided path exists and the amber.sh file exists
        if [ -f "${amber_path}/amber.sh" ]; then
            # Source the Amber environment
            source "${amber_path}/amber.sh"
            
            # Check if sourcing was successful
            if [ $? -ne 0 ]; then
                echo "Failed to source Amber. Please check the path."
                continue
            fi
            echo "Amber environment sourced successfully."
            break
        else
            echo "The provided path is incorrect or amber.sh is missing. Please try again."
        fi
    fi
done

# Function to get user input
get_user_input() {
    # Get total charge
    while true; do
        echo " "
        read -p "Please provide the total charge: " charge_total
        if [[ "$charge_total" =~ ^-?[0-9]+$ ]]; then
            break
        else
            echo "Invalid input. Please enter an integer for the total charge."
        fi
    done
    
    # Get input format for output of gaussian
    echo " "
    echo "Please select the input format:"
    echo "1- Gaussian output (log)"
    echo "2- Resp (gesp)"
    echo "3- PDB structure (pdb)"
    input_format=$(get_valid_input "Enter the number corresponding to the input format: " "1 2 3")

    # Map input format to appropriate flag
    case $input_format in
        1) input_form="gout";;
        2) input_form="gesp";;
        3) input_form="pdb";;
    esac

    # Get charge method
    echo " "
    echo "Please select the charge method (recommended: RESP): "
    echo "1- RESP (resp)"
    echo "2- Mulliken (mul)"
    echo "3- ESP (esp)"
    echo "4- AM1-BCC (bcc)"
    charge_method=$(get_valid_input "Enter the number corresponding to the charge method: " "1 2 3 4")

    # Map charge method to appropriate flag
    case $charge_method in
        1) method="resp";;
        2) method="mul";;
        3) method="esp";;
        4) method="bcc";;
    esac

    # Get atom type
    echo " "
    echo "Please select the atom type: "
    echo "1- Amber Force Field (AMBER)"
    echo "2- General Amber Force Field (GAFF)"
    echo "3- General Amber Force Field (GAFF2)"
    atom_type=$(get_valid_input "Enter the number corresponding to the atom type: " "1 2 3")

    # Map atom type to appropriate flag
    case $atom_type in
        1) at_type="amber";;
        2) at_type="gaff";;
        3) at_type="gaff2";;
    esac
}

# Call the function to get user input
get_user_input

# Function to execute antechamber with error handling
run_antechamber() {
    while true; do
	echo "      "
        read -p "Please provide the charge output file (e.g., .log, .resp, .pdb): " charge_data

        # Ensure the charge data file is found in RUN_DIR
        if [ ! -f "$RUN_DIR/$charge_data" ]; then
            echo "Charge data file not found in "$RUN_DIR". Please check the file name and try again."
            continue
        fi

        # Execute antechamber command with files from RUN_DIR
        antechamber -i "$RUN_DIR/$charge_data" -fi "$input_form" -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -c "$method" -s 2 -rn mol -nc "$charge_total" -j 5 -at "$at_type" -dr no > "$RUN_DIR/temp.dat"

        # Check if the command was successful
        if [ $? -ne 0 ]; then
            echo "Antechamber command failed. It seems there is no atomic charge calculation in the provided output file."
            retry=$(get_valid_input "Would you like to provide a different charge output file? (y/n): " "y n")

            case "$retry" in
                [yY])
                    echo "Retrying with a different charge output file..."
                    continue
                    ;;
                [nN])
                    echo "Exiting due to failed antechamber command."
                    exit 1
                    ;;
            esac
        else
            echo "Antechamber command executed successfully."
            break
        fi
    done
}

# Run the antechamber command
run_antechamber

# Ask the user for the checkpoint file
while true; do
    echo "      "
    read -p "Please provide the checkpoint file (.chk or .fchk): " checkpoint
    # Check if the file has a .chk or .fchk extension
    if [[ ! "$checkpoint" =~ \.(chk|fchk)$ ]]; then
        echo "Error: The file must have a .chk or .fchk extension. Please try again."
        continue
    fi
    # Check if the file exists in the RUN_DIR
    if [ ! -f "$RUN_DIR/$checkpoint" ]; then
        echo "Error: Checkpoint file not found in $RUN_DIR. Please check the file name and try again."
        continue
    fi

    # Get the file extension
    file_extension="${checkpoint##*.}"

    if [ "$file_extension" = "chk" ]; then
        while true; do
            echo "      "
            echo "Select your option:"
            echo "1- Gaussian is already loaded (formchk is available)"
            echo "2- Provide the Gaussian path"
            gauss_choice=$(get_valid_input "Enter your choice (1 or 2): " "1 2")
            formchk_command=""
            if [ "$gauss_choice" -eq 1 ]; then
                formchk_command="formchk"
            elif [ "$gauss_choice" -eq 2 ]; then
                echo "      "
                read -p "Please provide the path for Gaussian: " gauss_path
                formchk_path=$(find "$gauss_path" -type f -name "formchk" 2>/dev/null | head -n 1)
                if [ -z "$formchk_path" ]; then
                    echo "Error: Failed to find formchk in the provided path. Please check the path."
                    continue
                fi
                formchk_command="$formchk_path"
            else
                echo "Error: Invalid option selected."
                continue
            fi
            # Run formchk with the checkpoint file and output to complex.fchk
            if ! $formchk_command "$RUN_DIR/$checkpoint" "$RUN_DIR/complex.fchk" > "$RUN_DIR/temp.dat" 2>&1; then
                echo "Error: formchk command failed. Please check the input and try again."
                continue
            fi
            break
        done
    elif [ "$file_extension" = "fchk" ]; then
        # If the file is already .fchk, just copy it to complex.fchk
        cp "$RUN_DIR/$checkpoint" "$RUN_DIR/tempz.fchk"
        cp "$RUN_DIR/tempz.fchk" "$RUN_DIR/complex.fchk"
    fi

    break
done

# Ask the user to provide the optimized XYZ geometry file
while true; do
    echo "      "
    read -p "Please provide the optimized XYZ geometry file: " xyz_file
    
    # Check if the file has a .xyz extension
    if [[ ! "$xyz_file" =~ \.xyz$ ]]; then
        echo "Error: The file must have a .xyz extension. Please try again."
        continue
    fi
    
    # Check if the file exists in the RUN_DIR
    if [ ! -f "$RUN_DIR/$xyz_file" ]; then
        echo "Error: Optimized structure is not found in $RUN_DIR. Please check the file name and try again."
        continue
    fi
    
    break
done

# Executw the Python script to calculate bond angles and dihedrals
if [ -f "$SCRIPT_DIR/02_get_bond_angle.py" ]; then
    python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/$xyz_file"
else
    echo "Script 02_get_bond_angle.py not found in $SCRIPT_DIR. Exiting."
    exit 1
fi

# Check if the script was successful
if [ $? -ne 0 ]; then
    echo "Failed to execute 02_get_bond_angle.py. Exiting."
    exit 1
fi

# Run the Python and shell scripts in sequence with checks
# Generate a .mol2 file
# Correct the .mol2 file
# Generate a .frcmod file
# Assign new atom types in the .mol2 and .frcmod file
for script in 03_correct_mol2.py 04_parmch2_frcmod.sh 05_prepare_mol2_frcmod.py; do
    if [[ $script == *.py ]]; then
        if [ -f "$SCRIPT_DIR/$script" ]; then
            python3 "$SCRIPT_DIR/$script" "$RUN_DIR"
            if [ $? -ne 0 ]; then
                echo "Failed to execute $script. Exiting."
                exit 1
            fi
        else
            echo "Script $script not found in $SCRIPT_DIR. Exiting."
            exit 1
        fi
    elif [[ $script == *.sh ]]; then
        if [ -f "$SCRIPT_DIR/$script" ]; then
            bash "$SCRIPT_DIR/$script" "$RUN_DIR"
            if [ $? -ne 0 ]; then
                echo "Failed to execute $script. Exiting."
                exit 1
            fi
        else
            echo "Script $script not found in $SCRIPT_DIR. Exiting."
            exit 1
        fi
    fi
done

# Execute the Seminario method to calculate bond, angle, and force constant parameters.
python3 "$SCRIPT_DIR/Seminario_method.py" "$RUN_DIR/complex.fchk" > "$RUN_DIR/temp.dat"

if [ $? -ne 0 ]; then
    echo "Failed to execute Seminario_method.py. Exiting."
    exit 1
fi

# Generate the final paramter with the correct information
for script in 06_get_atom_type.py 07_Seminario_forcefield.py 08_update_forcefield.py 09_clean_updatedforcefield.py 10_postclean_updatedforcefield.py 11_retrieve_uffdata.py 13_final_clean.py; do
    if [ -f "$SCRIPT_DIR/$script" ]; then
        python3 "$SCRIPT_DIR/$script" "$RUN_DIR"
        if [ $? -ne 0 ]; then
            echo "Failed to execute $script. Exiting."
            exit 1
        fi
    else
        echo "Script $script not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi
done

# Change the name of outpur
mv filtered_COMPLEX_modified2.frcmod COMPLEX.frcmod
mv NEW_COMPLEX.mol2 COMPLEX.mol2 

# Prepare library input for tleap
echo "source leaprc.gaff" > input_library.tleap
echo "mol = loadmol2 "COMPLEX.mol2"" >> input_library.tleap
echo "loadamberparams COMPLEX.frcmod" >> input_library.tleap
echo "check mol" >> input_library.tleap
echo "charge mol" >> input_library.tleap
echo "savepdb mol COMPLEX.pdb" >> input_library.tleap
echo "saveoff mol COMPLEX.lib" >> input_library.tleap
echo "quit" >> input_library.tleap

tleap -f input_library.tleap > leap.log

#Function to correct lib file

python3 "$SCRIPT_DIR/12_generate_lib.py"

# Remove the unnecessary file
rm dihedral.dat distance.dat esout atom_type.dat COMPLEX_modified.mol2 COMPLEX_modified.frcmod complex.fchk forcefield2.dat forcefield.dat metal_number.dat new_atomtype.dat temp_COMPLEX_modified.frcmod temp.dat updated_COMPLEX_modified.frcmod updated_COMPLEX_modified2.frcmod angle.dat bond_angle_dihedral_data.dat qout punch QOUT ANTECHAMBER* ATOMTYPE.INF leap.log updated_updated_COMPLEX_modified2.frcmod 

# Print all the output name
echo ""
echo "Here is your mol2 COMPLEX.mol2 "
echo "Here is your frcmod COMPLEX.frcmod "
echo "Here is your pdb COMPLEX.pdb "
echo "Here is your lib COMPLEX.lib "

# Ask if the user wants to restrain the charge on specific atoms
echo " "
read -p "Would you like to restrain the charge on specific atoms? (Yes or No): " restrain_choice
if [[ "$restrain_choice" =~ ^[Yy][Ee][Ss]|[Yy]$ ]]; then
    echo " "
    while true; do
        num_atoms=$(get_valid_input "How many atoms do you want to restrain? " "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19")
        
        if [[ "$num_atoms" -lt 1 ]]; then
            echo "Error: The number of atoms to restrain must be at least 1. Please try again."
            continue
        fi

        > "$RUN_DIR/restraint_info.dat"  # Initialize the file to be empty

        for (( i=1; i<=num_atoms; i++ ))
        do
            while true; do
                echo " "
                read -p "Please provide the atom number and its charge for atom $i (e.g., 12 -0.834): " atom_number charge

                # Validate input
                if [[ ! "$atom_number" =~ ^[0-9]+$ ]] || [[ ! "$charge" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
                    echo "Error: Invalid input. Atom number must be an integer and charge must be a number."
                    continue
                fi

                echo "$atom_number $charge" >> "$RUN_DIR/restraint_info.dat"
                break
            done
        done

        # Ask if the user wants to review or redo the restraints
	echo " "
        read -p "Do you want to review the restraints? (y/n): " review_choice
        if [[ "$review_choice" =~ ^[Yy]$ ]]; then
            cat "$RUN_DIR/restraint_info.dat"
	    echo " "
            read -p "Do you want to redo the restraints? (y/n): " redo_choice
            if [[ "$redo_choice" =~ ^[Yy]$ ]]; then
                continue
            fi
        fi

        break
    done
    # Ask the user for the charge output
    while true; do
	echo " "
        read -p "Please enter the name of Gaussian output that contains all the ESP charges: " ESP_FILE
        
        # Check if the file has an accepted extension
        if [[ ! "$ESP_FILE" =~ \.(dat|log|out|txt)$ ]]; then
            echo "Error: The file must have a .log or .out extension. Please try again."
            continue
        fi
        
        # Check if the ESP file exists
        if [ ! -f "$RUN_DIR/$ESP_FILE" ]; then
            echo "Error: The Gaussian file '$ESP_FILE' was not found in $RUN_DIR. Please check the file name and try again."
            continue
        fi
        
        break  # Exit the loop if the file exists and has a valid extension
    done

    file_name="resp.dat"

    if [ -f "$file_name" ]; then
    	rm "$file_name"
    fi 

    # Extract relevant data from the Gaussian output
    grep "Atomic Center " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/a"
    grep "ESP Fit" "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/b"
    grep "Fit    " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/c"

    # detect the number of atoms
    i=$(wc -l < "$RUN_DIR/a")
    
    # Ensure the readit binary is in the script directory and is executable
    if [ ! -x "$SCRIPT_DIR/readit" ]; then
        echo "Error: readit binary not found or not executable in $SCRIPT_DIR"
        exit 1
    fi
    
    # Run the readit binary
    (cd "$RUN_DIR" && "$SCRIPT_DIR/readit")
    
    # Generate the required input for resp
    echo "floBF-resp run #1" > complex_esp.in
    echo " &cntrl" >> complex_esp.in 
    echo " nmol=1," >> complex_esp.in
    echo " ihfree=1," >> complex_esp.in
    echo " qwt=0.0005," >> complex_esp.in
    echo " iqopt=2," >> complex_esp.in
    echo " /" >> complex_esp.in
    echo "    1.0" >> complex_esp.in
    echo "Complex" >> complex_esp.in
    echo "   $charge_total   $i" >> complex_esp.in

    OUTPUT_FILE=complex_esp.in

# Read the atom numbers and charges from restrain into
    declare -A restrained_atoms
    while read -r atom_number charge; do
    	restrained_atoms[$atom_number]=$charge
    done < "$RUN_DIR/restraint_info.dat"

# Read the similar.dat file to organize atoms with similar chemical environments
    declare -A atom_types
    while read -r atom_number atom_type; do
    	atom_types[$atom_number]=$atom_type
    done < "$RUN_DIR/similar.dat"

# Convert arrays to space-separated lists
    restrained_atoms_list=$(for key in "${!restrained_atoms[@]}"; do echo "$key ${restrained_atoms[$key]}"; done)
    atom_types_list=$(for key in "${!atom_types[@]}"; do echo "$key ${atom_types[$key]}"; done)

# Extract and process the data relsted to atomic number, total charge, atom id
    awk -v i="$i" -v OFILE="$OUTPUT_FILE" -v restrained_atoms_list="$restrained_atoms_list" -v atom_types_list="$atom_types_list" '
BEGIN {
    start = 0
    # Populate associative arrays in awk
    split(restrained_atoms_list, restrained_array, " ")
    for (i in restrained_array) {
        restrained_atoms[restrained_array[i]] = 1
    }
    split(atom_types_list, atom_types_array, " ")
    for (i = 1; i <= length(atom_types_array); i+=2) {
        atom_types[atom_types_array[i]] = atom_types_array[i+1]
    }
}
/Standard orientation:/ { 
    start = 1
    getline; getline; getline; getline  # Skip 4 lines after "Standard orientation:"
    next
}
start && NF == 6 && $1 ~ /^[0-9]+$/ {
    atom_num = $1
    atom_type = $3
    if (atom_num in restrained_atoms) {
        atom_type = -1
    } else if (atom_num in atom_types) {
        atom_type = atom_types[atom_num]
    }
    printf "%3d%4d\n", $2, atom_type >> OFILE
}
/^---------------------------------------------------------------------$/ && start {
    exit
}
' "$RUN_DIR/$ESP_FILE"
    
    echo " " >> complex_esp.in
    echo " " >> complex_esp.in
    # Fill the input with the correct info about the charge
    awk -v lines="$i" 'BEGIN {
        cols = 8;  # Number of columns per line
        for (n = 1; n <= lines; n++) {
            printf "0.000000";
            if (n % cols == 0) {
                print "";  # New line after every 8th column
            } else if (n < lines) {
                printf "  ";  # Two spaces between charges
            }
        }
        if (lines % cols != 0) {
            print "";  # New line if last line has less than 8 columns
        }
    }' > "$RUN_DIR/complex_esp.qin"

    # Update complex_esp.qin by replacing the correct positions with their corresponding charges
    awk -v input="$RUN_DIR/restraint_info.dat" -v cols=8 -v total="$i" '
    BEGIN {
        # Initialize charges array with default zeros
        for (n = 1; n <= total; n++) {
            charges[n] = "0.000000";
        }

        # Replace the values at specific positions with the charges from restraint_info.dat
        while ((getline < input) > 0) {
            atom_number = $1;
            charge = sprintf("%.6f", $2);
            if (atom_number > 0 && atom_number <= total) {
                charges[atom_number] = charge;
            }
        }

        # Print the charges distributed across multiple lines with 8 columns
        for (n = 1; n <= total; n++) {
            printf "%s", charges[n];
            if (n % cols == 0) {
                print "";  # New line after every 8th column
            } else if (n < total) {
                printf "  ";  # Two spaces between charges
            }
        }
        if (total % cols != 0) {
            print "";  # New line if last line has less than 8 columns
        }
    }' > "$RUN_DIR/tmpfile" && mv "$RUN_DIR/tmpfile" "$RUN_DIR/complex_esp.qin"
    # Excute the resp from ambertools
    resp -O -i complex_esp.in -o complex_esp.out -p complex_esp.pch -t complex_esp.chg -q complex_esp.qin -e esp.dat

    # Reshape the charges in complex_esp.chg to be one column and save to charges.dat
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/COMPLEX_charged.mol2"

    awk '{for (i=1; i<=NF; i++) print $i}' "$RUN_DIR/complex_esp.chg" > "$RUN_DIR/charges.dat"

    # Replace the charges in COMPLEX_charged.mol2 with those from charges.dat
    awk -v charge_file="$RUN_DIR/charges.dat" '
    BEGIN {
        # Read all charges from charges.dat into an array
        i = 1;
        while ((getline line < charge_file) > 0) {
            charges[i] = line;
            i++;
        }
    }
    # Identify the section of the mol2 file containing the atoms
    /@<TRIPOS>ATOM/ {
        in_atoms = 1;
        atom_index = 1;
        print;
        next;
    }
    /@<TRIPOS>BOND/ {
        in_atoms = 0;
        print;
        next;
    }
    in_atoms {
        # Replace the 9th column with the corresponding charge from charges.dat
        $9 = charges[atom_index];
        atom_index++;
        # Print the line with fixed formatting to maintain original alignment
        printf "%7d %-8s %8.4f %8.4f %8.4f %-6s %4d %-5s %8.6f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9;
        next;
    }
    {
        print;  # Print other lines unchanged
    }
    ' "$RUN_DIR/COMPLEX_charged.mol2" > "$RUN_DIR/tmpfile" && mv "$RUN_DIR/tmpfile" "$RUN_DIR/COMPLEX_charged.mol2"
    
    #Generate copies for the files
    cp COMPLEX.mol2 oldCOMPLEX.mol2
    cp COMPLEX.lib  oldCOMPLEX.lib

    cp COMPLEX_charged.mol2 COMPLEX.mol2
   
    tleap -f input_library.tleap > leap.log

    #Function to correct lib file
    
    python3 "$SCRIPT_DIR/12_generate_lib.py"
    
    cp COMPLEX.lib COMPLEX_charged.lib
    cp oldCOMPLEX.mol2 COMPLEX.mol2
    cp oldCOMPLEX.lib COMPLEX.lib
    # Print output name
    echo " "
    echo "Charges updated in COMPLEX_charged.mol2"
    echo "lib file updated in COMPLEX_charged.lib"

    # Clean up temporary files
    rm -f "$RUN_DIR/a" "$RUN_DIR/b" "$RUN_DIR/c" "$RUN_DIR/tmpfile" complex_esp.in complex_esp.out complex_esp.qin complex_esp.chg complex_esp.pch esp.dat charges.dat restraint_info.dat oldCOMPLEX.mol2 oldCOMPLEX.lib leap.log
else
    echo "No charge restraints applied. Exiting."
fi

rm similar.dat input_library.tleap distance_type.dat tempz.fchk 
