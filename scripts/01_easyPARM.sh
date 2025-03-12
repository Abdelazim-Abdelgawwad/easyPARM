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


# Detect the directory where the current script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Get the directory from which easyPARM.sh was executed 
RUN_DIR="$1"

echo " "
echo "Running in directory: $RUN_DIR"
echo "Script directory: $SCRIPT_DIR"

# Ensure all outputs go to RUN_DIR
cd "$RUN_DIR"
# Remove files if it exists
files_to_remove=("dihedral.dat" "distance.dat" "esout" "atom_type.dat" "COMPLEX_modified.mol2" "COMPLEX_modified.frcmod" \
    "forcefield2.dat" "forcefield.dat" "metal_number.dat" "new_atomtype.dat" "temp_COMPLEX_modified.frcmod" \
    "temp.dat" "updated_COMPLEX_modified.frcmod" "updated_COMPLEX_modified2.frcmod" "angle.dat" "bond_angle_dihedral_data.dat" "new_atomtype1.dat"\
    "qout" "punch" "QOUT" "ATOMTYPE.INF" "leap.log" "updated_updated_COMPLEX_modified2.frcmod" "metals_complete.dat" "more_metal.dat" "new_atomtype2.dat" "REF_COMPLEX.mol2" "limited_data.dat" "line_number.dat" "ONE.mol2" "Reference_atom_type.dat" "COMPLEX.lib") 

for file in "${files_to_remove[@]}"; do
    if [ -e "$file" ]; then
        rm -f "$file"
    fi
done

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
echo " "
echo "================================="
echo "  Amber Configuration Menu "
echo "================================="
echo "Select your option:"
echo "1- Use currently loaded Amber"
echo "2- Specify Amber installation path"

while true; do
    choice=$(get_valid_input "Enter your choice: " "1 2")
    
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
	    echo " "
            echo "Invalid input. Please enter an integer for the total charge."
        fi
    done

    while true; do
	echo " "
        read -p "Please provide the total multiplicity: " multi_total
        if [[ "$multi_total" =~ ^-?[0-9]+$ ]]; then
            break
        else
	    echo " "
            echo "Invalid input. Please enter an integer for the total multiplicity."
        fi
    done
    # Ask the user to provide the optimized XYZ geometry file
    while true; do
	echo " "
        read -p "Please provide the optimized XYZ geometry file: " xyz_file
        if [[ ! "$xyz_file" =~ \.xyz$ ]]; then
            echo "Error: The file must have a .xyz extension. Please try again."
            continue
        fi
        if [ ! -f "$RUN_DIR/$xyz_file" ]; then
            echo "Error: Optimized structure is not found in $RUN_DIR. Please check the file name and try again."
            continue
        fi
        break
    done

    # Execute the Python script to calculate bond, angles and dihedrals
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

    echo " "
    echo "================================="
    echo "  Charge Method Selection Menu "
    echo "================================="
    echo "Select the charge calculation method:"
    echo "1- GAUSSIAN (RESP Charges)"
    echo "2- ORCA (CHELPG Charges)"
    echo "3- ORCA (RESP Charges)"
    echo "4- GAMESS (RESP Charges)"
    echo "5- GAMESS (GAMESS Fit Charges)"
    charge_output=$(get_valid_input "Enter your choice: " "1 2 3 4 5")

    if [ "$charge_output" -eq 1 ]; then
        get_gaussian_input
    elif [ "$charge_output" -eq 2 ]; then
        get_orca_input
    elif [ "$charge_output" -eq 3 ]; then
        get_orca_input
    elif [ "$charge_output" -eq 4 ]; then
        get_gamess_input
    elif [ "$charge_output" -eq 5 ]; then
        get_gamess_input
    fi
}

get_gaussian_input() {
    # Get input format for output of gaussian
    echo " "
    echo "================================="
    echo "  Charge Output Selection Menu "
    echo "================================="
    echo "Please select the input format:"
    echo "1- Gaussian Output (.log file)"
    echo "2- Resp (.gesp file)"
    input_format=$(get_valid_input "Enter your choice: " "1 2")

    # Map input format to appropriate flag
    case $input_format in
        1) input_form="gout";;
        2) input_form="gesp";;
    esac

    # Get charge method
    echo " "
    echo "================================="
    echo "  Charge Method Selection Menu "
    echo "================================="
    echo "Please select the charge method (recommended: RESP): "
    echo "1- RESP (resp)"
    echo "2- Mulliken (mul)"
    echo "3- ESP (esp)"
    echo "4- AM1-BCC (bcc)"
    charge_method=$(get_valid_input "Enter your choice: " "1 2 3 4")

    # Map charge method to appropriate flag
    case $charge_method in
        1) method="resp";;
        2) method="mul";;
        3) method="esp";;
        4) method="bcc";;
    esac

    get_atom_type
    run_antechamber_gaussian
}

get_orca_input() {
    get_atom_type
    run_antechamber_orca
}

get_gamess_input() {
    get_atom_type
    run_antechamber_gamess
}


get_atom_type() {
    echo " "
    echo "================================="
    echo "  Atom Type Selection Menu "
    echo "================================="
    echo "Please select the atom type: "
    echo "1- Amber Force Field (AMBER)"
    echo "2- General Amber Force Field (GAFF)"
    echo "3- General Amber Force Field (GAFF2)"
    atom_type=$(get_valid_input "Enter your choice: " "1 2 3")

    # Map atom type to appropriate flag
    case $atom_type in
        1) at_type="amber";;
        2) at_type="gaff";;
        3) at_type="gaff2";;
    esac
}

run_antechamber_gaussian() {
    while true; do
        echo " "
        read -p "Please provide the charge output file (e.g., .log, .gesp): " charge_data
        
        # Ensure the charge data file is found in RUN_DIR
        if [ ! -f "$RUN_DIR/$charge_data" ]; then
            echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
            continue
        fi

        # First antechamber attempt
        antechamber -i "$RUN_DIR/$charge_data" -fi "$input_form" -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -c "$method" -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -j 5 -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1 
	sed -i 's/ c2 / cc /g' $RUN_DIR/COMPLEX.mol2
	sed -i 's/ n2 / nd /g' $RUN_DIR/COMPLEX.mol2
        # Convert XYZ to PDB
        python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/$xyz_file" "$RUN_DIR/COMPLEX.pdb"
        
        # Revise Atom Type
        python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"

        # Additional processing steps to ensure that COMPLEX.mol2 is correct
        if [ ! -f "$RUN_DIR/limited_data.dat" ]; then
            if [ -f "$RUN_DIR/no_metal.dat" ]; then
                antechamber -i "$RUN_DIR/$charge_data" -fi "$input_form" -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -c "$method" -s 2 -rn mol -nc "$charge_total" -at "$at_type" -m "$multi_total" > "$RUN_DIR/temp.dat" 2>&1
            else
                python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
                antechamber -i "$RUN_DIR/mol.pdb" -fi pdb -o "$RUN_DIR/ONE.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
                python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
                if [ -f "$RUN_DIR/ONE.mol2" ]; then
                    mv "$RUN_DIR/COMPLEX_modified.mol2" "$RUN_DIR/COMPLEX.mol2"
                fi
            fi
        else 
            break 
        fi

        # Check if COMPLEX.mol2 was generated after all attempts
        if [ ! -f "$RUN_DIR/COMPLEX.mol2" ]; then
            echo "  "
	    echo "Failed to generate COMPLEX.mol2 after multiple attempts."
            retry=$(get_valid_input "Would you like to provide a different charge output file? (y/n)" "y n")
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
            echo "Antechamber command executed successfully. "
            break
        fi
    done
}

run_antechamber_orca() {
    # Convert xyz to pdb
    if [ -f "$SCRIPT_DIR/xyz_to_pdb.py" ]; then
        python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/$xyz_file" "$RUN_DIR/COMPLEX.pdb"
    else
        echo "Script xyz_to_pdb.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi
    
    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute xyz_to_pdb.py. Exiting."
        exit 1
    fi
    
    # Read charge output
    while true; do
        echo " "
        if [ "$charge_output" -eq 2 ]; then
		read -p "Please provide the charge output file (e.g., .log, .out): " charge_data
		if [[ ! "$charge_data" =~ \.(log|out)$ ]]; then
			echo "Error: The file must have a .log or .out extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi
	elif [ "$charge_output" -eq 3 ]; then
		read -p "Please provide the charge output file ( .vpot): " charge_data
		if [[ ! "$charge_data" =~ \.vpot$ ]]; then
			echo "Error: The file must have a .vpot extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi

	fi
        break
    done
    
    # First antechamber attempt
    antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no -j 5 > "$RUN_DIR/temp.dat" 2>&1
    sed -i 's/ c2 / cc /g' $RUN_DIR/COMPLEX.mol2
    sed -i 's/ n2 / nd /g' $RUN_DIR/COMPLEX.mol2
    
    # Revise Atom Type
    python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
    
    # Additional processing steps
    if [ ! -f "$RUN_DIR/limited_data.dat" ]; then
        if [ -f "$RUN_DIR/no_metal.dat" ]; then
            antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" > "$RUN_DIR/temp.dat" 2>&1
        else 
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            antechamber -i "$RUN_DIR/mol.pdb" -fi pdb -o "$RUN_DIR/ONE.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            if [ -f "$RUN_DIR/ONE.mol2" ]; then
                mv "$RUN_DIR/COMPLEX_modified.mol2" "$RUN_DIR/COMPLEX.mol2"
            fi
        fi
    else 
        break 
    fi
    
    # Check if COMPLEX.mol2 was generated after all attempts
    if [ ! -f "$RUN_DIR/COMPLEX.mol2" ]; then
        echo " "
	echo "Failed to generate COMPLEX.mol2 after multiple attempts."
        retry=$(get_valid_input "Would you like to provide a different charge output file? (y/n)" "y n")
        case "$retry" in
            [yY]) 
                echo "Retrying with a different charge output file..."
                run_antechamber_orca
                return
                ;;
            [nN]) 
                echo "Exiting due to failed antechamber command."
                exit 1
                ;;
        esac
    else
        echo "Antechamber command executed successfully. "
    fi

    # Update mol2 file with CHELPG charges (if selected)
    if [ "$charge_output" -eq 2 ]; then
        awk '/CHELPG Charges/,/Total charge:/' "$RUN_DIR/$charge_data" | \
            grep -E '^\s*[0-9]+' | \
            awk '{print $NF}' > "$RUN_DIR/charges.dat"
        cp "$RUN_DIR/charges.dat" "$RUN_DIR/charges.chg"
       	python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
    
    fi
}
run_antechamber_gamess() {
    # Convert xyz to pdb
    if [ -f "$SCRIPT_DIR/xyz_to_pdb.py" ]; then
        python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/$xyz_file" "$RUN_DIR/COMPLEX.pdb"
    else
        echo "Script xyz_to_pdb.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi
    
    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute xyz_to_pdb.py. Exiting."
        exit 1
    fi
    
    # Read charge output
    while true; do
	echo " "
        if [ "$charge_output" -eq 4 ]; then
		read -p "Please provide the charge output file ( .dat): " charge_data
		if [[ ! "$charge_data" =~ \.dat$ ]]; then
			echo "Error: The file must have a .dat extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi
	elif [ "$charge_output" -eq 5 ]; then
		read -p "Please provide the charge output file ( .log): " charge_data
		if [[ ! "$charge_data" =~ \.log$ ]]; then
			echo "Error: The file must have a .log extension. Please try again."
			continue
		fi
		if [ ! -f "$RUN_DIR/$charge_data" ]; then
			echo "Charge data file not found in $RUN_DIR. Please check the file name and try again."
			continue
		fi

	fi
        break
    done
 
    # First antechamber attempt
    antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -at "$at_type" -m "$multi_total" -dr no -j 5 > "$RUN_DIR/temp.dat" 2>&1
    sed -i 's/ c2 / cc /g' $RUN_DIR/COMPLEX.mol2
    sed -i 's/ n2 / nd /g' $RUN_DIR/COMPLEX.mol2
      
    # Revise Atom Type
    python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
    
    # Additional processing steps
    if [ ! -f "$RUN_DIR/limited_data.dat" ]; then
        if [ -f "$RUN_DIR/no_metal.dat" ]; then
            antechamber -i "$RUN_DIR/COMPLEX.pdb" -fi pdb -o "$RUN_DIR/COMPLEX.mol2" -fo mol2 -s 2 -rn mol -m "$multi_total" -nc "$charge_total" -at "$at_type" > "$RUN_DIR/temp.dat" 2>&1
        else 
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            antechamber -i "$RUN_DIR/mol.pdb" -fi pdb -o "$RUN_DIR/ONE.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
            python3 "$SCRIPT_DIR/Revise_Atom_Type.py" > "$RUN_DIR/temp.dat"
            if [ -f "$RUN_DIR/ONE.mol2" ]; then
                mv "$RUN_DIR/COMPLEX_modified.mol2" "$RUN_DIR/COMPLEX.mol2"
            fi
        fi
    else 
        break 
    fi
    
    # Check if COMPLEX.mol2 was generated after all attempts
    if [ ! -f "$RUN_DIR/COMPLEX.mol2" ]; then
        echo " "
	echo "Failed to generate COMPLEX.mol2 after multiple attempts."
        retry=$(get_valid_input "Would you like to provide a different charge output file? (y/n)" "y n")
        case "$retry" in
            [yY]) 
                echo "Retrying with a different charge output file..."
                run_antechamber_orca
                return
                ;;
            [nN]) 
                echo "Exiting due to failed antechamber command."
                exit 1
                ;;
        esac
    else
        echo "Antechamber command executed successfully. "
    fi
    
    # Update mol2 file with gamess charges (if selected)
    if [ "$charge_output" -eq 5 ]; then
	awk '/NET CHARGES:/,/RMS DEVIATION IS/ {
		if ($1 ~ /^[A-Z]/ && $3 ~ /^-?[0-9]/)
    			print $2
	}' "$RUN_DIR/$charge_data" > "$RUN_DIR/charges.dat"
        cp "$RUN_DIR/charges.dat" "$RUN_DIR/charges.chg"
       	python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
    fi

}

# Main execution
get_user_input

# Ask the user for the input file
echo " "
echo "================================="
echo "  Input Format Selection Menu "
echo "================================="
echo "Please select the format you will provide:"
echo "1- Orca Output"
echo "2- Gaussian Output"
echo "3- Gaussian Checkpoint"
echo "4- Gaussian Formatted Checkpoint"
echo "5- Gamess Output"

qm_output=$(get_valid_input "Enter your choice: " "1 2 3 4 5")

case $qm_output in
    1)
        echo " "
	echo "You've selected Orca Output"
	if [ "$charge_output" -eq 1 -o "$charge_output" -eq 3 -o "$charge_output" -eq 4 -o "$charge_output" -eq 5 ]; then
            while true; do
                echo " "
                read -p "Please provide the Orca output file (.log or .out): " orca_output
                if [ ! -f "$RUN_DIR/$orca_output" ]; then
                    echo "Orca output file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            while true; do
                echo " "
                read -p "Please provide the Orca hessian file (.hess): " orca_hessian
                if [ ! -f "$RUN_DIR/$orca_hessian" ]; then
                    echo "Orca hessian file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            python3 "$SCRIPT_DIR/Seminario_method_ORCA.py" "$RUN_DIR/$orca_hessian" "$RUN_DIR/$orca_output" > "$RUN_DIR/temp.dat"
        elif [ "$charge_output" -eq 2 ]; then
            while true; do
                echo " "
                read -p "Please provide the Orca hessian file (.hess): " orca_hessian
                if [ ! -f "$RUN_DIR/$orca_hessian" ]; then
                    echo "Orca hessian file not found in $RUN_DIR. Please check the file name and try again."
                    continue
                fi
                break
            done
            python3 "$SCRIPT_DIR/Seminario_method_ORCA.py" "$RUN_DIR/$orca_hessian" "$RUN_DIR/$charge_data" > "$RUN_DIR/temp.dat"
        fi
        ;;
    2)
        echo " "
        echo "You've selected Gaussian output"
        while true; do
            echo " "
            read -p "Please provide the Gaussian output file (.log or .out): " gaussian_output
            if [ ! -f "$RUN_DIR/$gaussian_output" ]; then
                echo "Gaussian output file not found in $RUN_DIR. Please check the file name and try again."
                continue
            fi
            break
        done
        # Execute the Seminario method to calculate bond, angle, and force constant parameters.
        python3 "$SCRIPT_DIR/Seminario_method_GAUSSIAN.py" "$RUN_DIR/$gaussian_output" 2 > "$RUN_DIR/temp.dat"
        if [ $? -ne 0 ]; then
            echo "Failed to execute Seminario_method_GAUSSIAN.py. Exiting."
            exit 1
        fi
        ;;
    3)
        echo " "
        echo "You've selected Gaussian checkpoint"
        while true; do
            echo " "
            read -p "Please provide the checkpoint file (.chk): " checkpoint
            if [[ ! "$checkpoint" =~ \.chk$ ]]; then
                echo "Error: The file must have a .chk extension. Please try again."
                continue
            fi
            if [ ! -f "$RUN_DIR/$checkpoint" ]; then
                echo "Error: Checkpoint file not found in $RUN_DIR. Please check the file name and try again."
                continue
            fi
            break
        done
        echo " "
        echo "================================="
        echo "  Gaussian Selection Menu "
        echo "================================="
        echo "1- Gaussian is already loaded (formchk is available)"
        echo "2- Provide the Gaussian path"
        gauss_choice=$(get_valid_input "Enter your choice: " "1 2")

        if [ "$gauss_choice" -eq 1 ]; then
            formchk_command="formchk"
        elif [ "$gauss_choice" -eq 2 ]; then
            while true; do
                read -p "Please provide the path for Gaussian: " gauss_path
                formchk_path=$(find "$gauss_path" -type f -name "formchk" 2>/dev/null | head -n 1)
                if [ -z "$formchk_path" ]; then
                    echo "Error: Failed to find formchk in the provided path. Please check the path and try again."
                    continue
                fi
                formchk_command="$formchk_path"
                break
            done
        fi

        # Run formchk
        if ! $formchk_command "$RUN_DIR/$checkpoint" "$RUN_DIR/complex.fchk" > "$RUN_DIR/temp.dat" 2>&1; then
            echo "Error: formchk command failed. Please check the input and try again."
            exit 1
        fi
        # Execute the Seminario method to calculate bond, angle, and force constant parameters.
        python3 "$SCRIPT_DIR/Seminario_method_GAUSSIAN.py" "$RUN_DIR/complex.fchk" 3 > "$RUN_DIR/temp.dat"
        if [ $? -ne 0 ]; then
            echo "Failed to execute Seminario_method_GAUSSIAN.py. Exiting."
            exit 1
        fi
        ;;
    4)
        echo " "
        echo "You've selected Gaussian formatted checkpoint"
        while true; do
            echo " "
            read -p "Please provide the formatted checkpoint file (.fchk): " fchk_file
            # Validation for .fchk file
            if [[ ! "$fchk_file" =~ \.fchk$ ]]; then
                echo "Error: The file must have a .fchk extension. Please try again."
                continue
            fi
            if [ ! -f "$RUN_DIR/$fchk_file" ]; then
                echo "Error: Formatted checkpoint file not found in $RUN_DIR. Please check the file name and try again."
                continue
            fi
            break
        done
        cp "$RUN_DIR/$fchk_file" "$RUN_DIR/temp.fchk"
        mv "$RUN_DIR/temp.fchk" "$RUN_DIR/complex.fchk"
        # Execute the Seminario method to calculate bond, angle, and force constant parameters.
        python3 "$SCRIPT_DIR/Seminario_method_GAUSSIAN.py" "$RUN_DIR/complex.fchk" 4 > "$RUN_DIR/temp.dat"
        if [ $? -ne 0 ]; then
            echo "Failed to execute Seminario_method_GAUSSIAN.py. Exiting."
            exit 1
        fi
        ;;
    5)
        echo " "
        echo "You've selected Gamess output"
        while true; do
            echo " "
            read -p "Please provide the Gamess output file ( .dat): " gamess_output
            if [ ! -f "$RUN_DIR/$gamess_output" ]; then
                echo "Gaussian output file not found in $RUN_DIR. Please check the file name and try again."
                continue
            fi
            break
        done
        # Execute the Seminario method to calculate bond, angle, and force constant parameters.
        python3 "$SCRIPT_DIR/Seminario_method_GAMESS.py" "$RUN_DIR/$gamess_output" > "$RUN_DIR/temp.dat"
        if [ $? -ne 0 ]; then
            echo "Failed to execute Seminario_method_GAMESS.py. Exiting."
            exit 1
        fi
        ;;
esac
    
if [ "$charge_output" -eq 3 ]; then
        # For ORCA vpot file, simply copy the file
        cp "$RUN_DIR/$charge_data" "$RUN_DIR/esp.in"
	# Read the first line, remove the space, and store it
        first_line=$(sed '1q;d' "$RUN_DIR/esp.in" | sed -E 's/^([[:space:]]+)([0-9]+)[[:space:]]+([0-9]+)/\1\2\3/')
	# Replace the first line in the file
        sed -i "1s/.*/$first_line/" "$RUN_DIR/esp.in"
        # Generate the required input for resp calculation  
	python3 "$SCRIPT_DIR/RESP_ORCA.py" "$RUN_DIR/$orca_hessian" $charge_total "$RUN_DIR/similar.dat" > "$RUN_DIR/temp.dat"
	# Run resp 
    	resp -O -i "$RUN_DIR/resp.in" -o "$RUN_DIR/resp.out" -e "$RUN_DIR/esp.in" -t "$RUN_DIR/esp.chg"
	# reshape charges
        awk '{for (i=1; i<=NF; i++) print $i}' "$RUN_DIR/esp.chg" > "$RUN_DIR/charges.chg"
	# update mol2 charges
        python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
elif [ "$charge_output" -eq 4 ]; then
        python3 "$SCRIPT_DIR/RESP_GAMESS.py" "$RUN_DIR/$charge_data" $charge_total "$RUN_DIR/similar.dat" > "$RUN_DIR/temp.dat"
    	resp -O -i "$RUN_DIR/resp.in" -o "$RUN_DIR/resp.out" -e "$RUN_DIR/esp.in" -t "$RUN_DIR/esp.chg"
	# reshape charges
        awk '{for (i=1; i<=NF; i++) print $i}' "$RUN_DIR/esp.chg" > "$RUN_DIR/charges.chg"
	# update mol2 charges
        python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
fi

#Remove unnecessary file
if [  -f "$RUN_DIR/metals_complete.dat" ]; then
    rm "$RUN_DIR/metals_complete.dat" 
fi

# Run the Python and shell scripts in sequence with checks
# Generate a .mol2 file
# Correct the .mol2 file
# Generate a .frcmod file
# Assign new atom types in the .mol2 and .frcmod file
for script in 03_correct_mol2.py 04_parmch2_frcmod.sh; do
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

echo " "
metalloprotein_choice=$(get_valid_input "Does your structure belong to MetalloProtein ? (y/n): " "y n yes no Y N YES NO Yes No")
if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
    echo " "
    read -p "Please provide the metalloprotein pdb file: " protein_pdb

    pdb4amber -i "$RUN_DIR/$protein_pdb" -o "$RUN_DIR/metalloprotein_easyPARM.pdb"  > "$RUN_DIR/temp.dat" 2>&1 
    
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/NEW_COMPLEX.mol2" 
    if [ -f "$SCRIPT_DIR/metalloprotein.py" ]; then
        python3 "$SCRIPT_DIR/metalloprotein.py" 
    	mv "$RUN_DIR/COMPLEX_updated.mol2" "$RUN_DIR/easyCOMPLEX.mol2"
    else
        echo "Script metalloprotein.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi

    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute metalloprotein.py. Exiting."
        exit 1
    fi
    mv "$RUN_DIR/easyCOMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
    python3 "$SCRIPT_DIR/distribute_metalloprotein_charge.py" "$RUN_DIR/$protein_pdb" "$RUN_DIR/COMPLEX_modified.mol2" $charge_total 
    mv "$RUN_DIR/updated_easy_COMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
else
    :
fi

# Select appropriate preparation script in order to generate the frcmod depend on the complexation of the system
# more_metal.dat refer that there are more than metal or non standard atom then the code will select the code depend on that
if [ -f "$RUN_DIR/more_metal.dat" ]; then

    python3 "$SCRIPT_DIR/05_prepare_mol2_frcmod_more_atom.py" "$RUN_DIR"
    cp "$RUN_DIR/new_atomtype.dat" "$RUN_DIR/new_atomtype1.dat"    
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/REF_COMPLEX.mol2"    
    cp "$RUN_DIR/NEW_COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"    
    mv "$RUN_DIR/metal_number.dat" "$RUN_DIR/metals_complete.dat"
    rm "$RUN_DIR/more_metal.dat"
    python3 "$SCRIPT_DIR/generate_preforcefield.py" "$RUN_DIR"
    if [ $? -ne 0 ]; then
        echo "Failed to execute 05_prepare_mol2_frcmod_more_atom.py. Exiting."
        exit 1
    fi
fi

if [ ! -f "$RUN_DIR/more_metal.dat" ]; then
    python3 "$SCRIPT_DIR/05_prepare_mol2_frcmod.py" "$RUN_DIR"
    if [ $? -ne 0 ]; then
        echo "Failed to execute 05_prepare_mol2_frcmod.py. Exiting."
        exit 1
    fi
fi

if [  -f "$RUN_DIR/metals_complete.dat" ]; then
    cp "$RUN_DIR/metals_complete.dat" "$RUN_DIR/metal_number.dat"
fi

if [  -f "$RUN_DIR/new_atomtype1.dat" ]; then
    cp "$RUN_DIR/new_atomtype.dat" "$RUN_DIR/new_atomtype2.dat"    
    cp "$RUN_DIR/REF_COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
    cat "$RUN_DIR/new_atomtype1.dat" >> "$RUN_DIR/new_atomtype2.dat"
fi


echo " "
if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
 
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/1COMPLEX.mol2"
    if [ -f "$SCRIPT_DIR/metalloprotein.py" ]; then
        python3 "$SCRIPT_DIR/metalloprotein.py" 
    	mv "$RUN_DIR/COMPLEX_updated.mol2" "$RUN_DIR/easyCOMPLEX.mol2"
    else
        echo "Script metalloprotein.py not found in $SCRIPT_DIR. Exiting."
        exit 1
    fi

    # Check if the script was successful
    if [ $? -ne 0 ]; then
        echo "Failed to execute metalloprotein.py. Exiting."
        exit 1
    fi
   
    python3 "$SCRIPT_DIR/xyz_to_pdb.py" "$RUN_DIR/part_QM.xyz" "$RUN_DIR/part_QM.pdb"
    antechamber -i "$RUN_DIR/part_QM.pdb" -fi pdb -o "$RUN_DIR/QM.mol2" -fo mol2 -s 2 -rn mol -nc "$charge_total" -m "$multi_total" -at "$at_type" -dr no > "$RUN_DIR/temp.dat" 2>&1
    sed -i '$d' "$RUN_DIR/nonstand.pdb"
    cat "$RUN_DIR/part_QM.pdb" >> "$RUN_DIR/nonstand.pdb"
    cat "$RUN_DIR/part_QM.pdb" >> "$RUN_DIR/easynonstands.pdb"
    pdb4amber -i "$RUN_DIR/easynonstands.pdb" -o "$RUN_DIR/easyPARM.pdb"  > "$RUN_DIR/temp.dat" 2>&1 
     
    input_file="$RUN_DIR/easyPARM_residues.dat"

    # Read each line from the input file
    while IFS=' ' read -r pdbout mol2out residue_name; do
    # Skip empty lines or lines starting with #
    	if [[ -z "$pdbout" || "$pdbout" == \#* ]]; then
        	continue
    	fi
    
    # Run the antechamber command
    	antechamber -i "$pdbout" -fi pdb -o "$mol2out" -fo mol2 -s 2 -rn "$residue_name" 
    done < "$input_file" > "$RUN_DIR/temp.dat" 2>&1

     
    python3 "$SCRIPT_DIR/update_metalloprotein_charge.py" 
    mv "$RUN_DIR/QM.mol2" "$RUN_DIR/COMPLEX.mol2"
    python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/part_QM.xyz"
    python3 "$SCRIPT_DIR/03_correct_mol2.py" 
    mv "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/QM.mol2"
    mv "$RUN_DIR/1COMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
    mv "$RUN_DIR/easyCOMPLEX.mol2" "$RUN_DIR/COMPLEX_modified.mol2" 
    python3 "$SCRIPT_DIR/02_get_bond_angle.py" "$RUN_DIR/$xyz_file"
else
    :
fi

# Generate the frcmod file with all the correct info and with a clean version
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

# Change the name of output
mv filtered_COMPLEX_modified2.frcmod COMPLEX.frcmod
cp NEW_COMPLEX.mol2 COMPLEX.mol2 

resid_ID=$(get_valid_input "Would you like to change the residue ID (Default= mol)? (y/n): " "y n yes no Y N YES NO Yes No")
    if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
	    read -r -p "Please provide the residue name: " resid_name
    fi

if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
    mv "$RUN_DIR/nonstand.pdb" "$RUN_DIR/easyPARM_MetalloProtein.pdb"
    mv "$RUN_DIR/QM.mol2" "$RUN_DIR/METAL.mol2"
    mv "$RUN_DIR/coordinated_residues.txt" "$RUN_DIR/Bond_Info.dat"
    awk '
    /^ANGLE/,/^DIHE/ {
        if ($2 ~ /[0-9]/) {
            if ($2 < 10) {
                $2 = $2 * 4.299
            } else if ($2 < 20) {
                $2 = $2 * 2.799
            } else if ($2 < 29) {
                $2 = $2 * 1.699
            }
        }
        print
        next
    }
    { print }' "$RUN_DIR/COMPLEX.frcmod" > "$RUN_DIR/temp.frcmod"
    mv "$RUN_DIR/temp.frcmod" "$RUN_DIR/COMPLEX.frcmod"
    echo "  "
    

    echo "=========================================="
    echo "         	Output Files "
    echo "=========================================="

    input_file="$RUN_DIR/easyPARM_residues.dat"
    line_number=1

# Read each line from the input file
    while IFS=' ' read -r pdbout mol2out residue_name; do
    	# Skip empty lines or lines starting with #
    	if [[ -z "$mol2out" || "$pdbout" == \#* ]]; then
        	continue
    	fi	
    
    # Print the mol2out with the line number
    	echo "Mol2 $line_number			  : $mol2out"
    	((line_number++))
    done < "$input_file" 
    echo "Mol2  		    	  : METAL.mol2"
    echo "Frcmod                    : COMPLEX.frcmod"
    echo "MetalloProtein pdb        : easyPARM_MetalloProtein.pdb"
    echo "Bond Information	  : Bond_Info.dat"
    echo "New Atom Type	  	  : Hybridization_Info.dat"
    if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
	    #sed -i -E "s/ mol /${resid_name}/g" "$RUN_DIR/COMPLEX.mol2"
	    sed -i "s/\<mol\>/${resid_name}/g" "$RUN_DIR/METAL.mol2"
    fi
else 
    	if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then
	    #sed -i -E "s/ mol /${resid_name}/g" "$RUN_DIR/COMPLEX.mol2"
	    	sed -i "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.mol2"
	    	sed -i "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.pdb"
    	fi

# Print all the output name
	echo "================================="
	echo "         Output Files "
	echo "================================="
	echo "Mol2   		: COMPLEX.mol2"
	echo "Frcmod 		: COMPLEX.frcmod"
	echo "PDB    		: COMPLEX.pdb"
	echo "Lib    		: COMPLEX.lib"
    	echo "New Atom Type	: Hybridization_Info.dat"

fi	

# Prepare library input for tleap
echo "source leaprc.gaff" > input_library.tleap
echo "loadamberparams COMPLEX.frcmod" >> input_library.tleap
echo "mol = loadmol2 "COMPLEX.mol2"" >> input_library.tleap
echo "check mol" >> input_library.tleap
echo "charge mol" >> input_library.tleap
echo "savepdb mol COMPLEX.pdb" >> input_library.tleap
echo "saveoff mol COMPLEX.lib" >> input_library.tleap
echo "quit" >> input_library.tleap

tleap -f input_library.tleap > leap.log

#Function to correct lib file

python3 "$SCRIPT_DIR/12_generate_lib.py"
if [[ "${resid_ID,,}" =~ ^(y|yes)$ ]]; then 
	sed -i "s/\<mol\>/${resid_name}/g" "$RUN_DIR/COMPLEX.lib"
fi

echo " "
charmm_FF=$(get_valid_input "Would you like to generate the CHARMM force field format? (y/n): " "y n yes no Y N YES NO Yes No")
if [[ "${charmm_FF,,}" =~ ^(y|yes)$ ]] && [[ "${metalloprotein_choice,,}" =~ ^(n|no)$ ]] ; then
	python3 "$SCRIPT_DIR/str_helper.py" > "$RUN_DIR/temp.dat" 2>&1 
	python3 "$SCRIPT_DIR/generate_str.py" "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/forcefield.dat" > "$RUN_DIR/COMPLEX.top"  
	python3 "$SCRIPT_DIR/assign_amber_charmm.py" > "$RUN_DIR/temp.dat" 2>&1 

	echo " "
	echo "================================="
	echo "      CHARMM Output Files "
	echo "================================="
	echo "Force Field 	: COMPLEX.prm"
	echo "PDB    		: COMPLEX.pdb"
	echo "Topology     	: COMPLEX.top"

fi


if [[ "${metalloprotein_choice,,}" =~ ^(y|yes)$ ]]; then
	rm COMPLEX.lib
	rm COMPLEX.mol2
        rm COMPLEX.pdb	
	rm metalloprotein_easyPARM_*
	rm *_renum.txt
	rm *_sslink
	rm *_nonprot.pdb
	rm charge_*.dat
	rm recalculated_charges.dat
	rm processed_charges.dat
	rm original_charges.dat
	rm coordination_analysis.txt
        rm charge_statistics.txt	
	rm fixed_charges.dat
fi
# Remove the unnecessary file
rm ANTECHAMBER*
files_to_remove=("dihedral.dat" "distance.dat" "esout" "atom_type.dat" "COMPLEX_modified.mol2" "COMPLEX_modified.frcmod" \
    "complex.fchk" "forcefield2.dat" "metal_number.dat" "temp_COMPLEX_modified.frcmod" "new_atomtype.dat" \
    "temp.dat" "updated_COMPLEX_modified.frcmod" "updated_COMPLEX_modified2.frcmod" "angle.dat" "new_atomtype1.dat"\
    "qout" "punch" "QOUT" "ATOMTYPE.INF" "leap.log" "updated_updated_COMPLEX_modified2.frcmod" "metals_complete.dat" "more_metal.dat" "new_atomtype2.dat" "REF_COMPLEX.mol2" "limited_data.dat" "mol.pdb" "line_number.dat" "ONE.mol2" "Reference_atom_type.dat" "REFQM.pdb" "NEW_COMPLEX.mol2"\
    "QM.pdb" "nonstand.pdb" "part_QM.xyz" "part_QM.pdb" "charge_qm.dat" "metalloprotein.pdb" "metalloprotein_easyPARM.pdb" "charges_all.dat" "easynonstands.pdb" "easyPARM.pdb" "easyPARM_residues.dat" ) 

for file in "${files_to_remove[@]}"; do
    if [ -e "$file" ]; then
        rm -f "$file"
    fi
done

# Ask if the user wants to restrain the charge on specific atoms
echo "  "
restrain_choice=$(get_valid_input "Would you like to restrain the charge on specific atoms? (y/n): " "y n yes no Y N YES NO Yes No")
if [[ "${restrain_choice,,}" =~ ^(y|yes)$ ]]; then
    echo " "
    echo "========================================================================"
    echo "                           Restrain Method "
    echo "========================================================================"
    while true; do
	echo "  "
        num_atoms=$(get_valid_input "How many atoms do you want to restrain? " "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 ")
        
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
        
	read -p "Please enter the name of output that contains all the ESP charges: " ESP_FILE
        
        # Check if the file has an accepted extension
        if [[ ! "$ESP_FILE" =~ \.(vpot|log|out|dat)$ ]]; then
            echo "Error: The file must have a .log or .dat or .vpot extension. Please try again."
            continue
        fi
        
        # Check if the ESP file exists
        if [ ! -f "$RUN_DIR/$ESP_FILE" ]; then
            echo "Error: The Gaussian file '$ESP_FILE' was not found in $RUN_DIR. Please check the file name and try again."
            continue
        fi
        
        break  # Exit the loop if the file exists and has a valid extension
    done

    file_name="esp.dat"

    if [ -f "$file_name" ]; then
    	rm "$file_name"
    fi 
    # Extract the file extension
extension="${ESP_FILE##*.}"

    # Check the extension and set qm_output
    if [ "$extension" = "log" ]; then
    	qm_output=2
    elif [ "$extension" = "vpot" ]; then
    	qm_output=1
    elif [ "$extension" = "dat" ]; then
    	qm_output=3
    else
    	echo "Error: Unrecognized file extension: .$extension"
    	exit 1
    fi
    
    if [ "$qm_output" -eq 1 ]; then
        # For ORCA log file, simply copy the file
        cp "$RUN_DIR/$ESP_FILE" "$RUN_DIR/esp.dat"
	# Read the first line, remove the space, and store it
        first_line=$(sed '1q;d' "$RUN_DIR/esp.dat" | sed -E 's/^([[:space:]]+)([0-9]+)[[:space:]]+([0-9]+)/\1\2\3/')
# Replace the first line in the file
        sed -i "1s/.*/$first_line/" "$RUN_DIR/esp.dat"
          
	i=$(wc -l < "$RUN_DIR/atomic_number.dat")

    elif [ "$qm_output" -eq 2 ] ; then
        # Extract relevant data from the output for gaussian
        grep "Atomic Center " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/a"
        grep "ESP Fit" "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/b"
        grep "Fit    " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/c"

        # Detect the number of atoms
        i=$(wc -l < "$RUN_DIR/a")

        # Ensure the readit binary is in the script directory and is executable
        if [ ! -x "$SCRIPT_DIR/readit" ]; then
            echo "Error: readit binary not found or not executable in $SCRIPT_DIR"
            exit 1
        fi

        # Run the readit binary
        (cd "$RUN_DIR" && "$SCRIPT_DIR/readit")

    elif [ "$qm_output" -eq 3 ] ; then
        # For GAMESS dat file, prepare resp files
        python3 "$SCRIPT_DIR/RESP_GAMESS.py" "$RUN_DIR/$ESP_FILE" $charge_total "$RUN_DIR/similar.dat" > "$RUN_DIR/temp.dat"
	cp "$RUN_DIR/esp.in" "$RUN_DIR/esp.dat"
	i=$(wc -l < "$RUN_DIR/atomic_number.dat")

    else
        exit 1
    fi

if [ "$qm_output" = "1" ] || [ "$qm_output" = "3" ]; then    
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
    # Process the atomic numbers and add them to complex_esp.in
    awk -v OFILE="$OUTPUT_FILE" -v restrained_atoms_list="$restrained_atoms_list" -v atom_types_list="$atom_types_list" '
    BEGIN {
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
    {
        atom_num = NR
        atomic_number = $1
        if (atom_num in restrained_atoms) {
            atom_type = -1
        } else if (atom_num in atom_types) {
            atom_type = atom_types[atom_num]
        } else {
            atom_type = 0
        }
        printf "%3d%4d\n", atomic_number, atom_type >> OFILE
    }
    ' atomic_number.dat
    echo " " >> complex_esp.in
    echo " " >> complex_esp.in

elif [ "$qm_output" = "2" ] ; then
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
    # Extract and process the data related to atomic number, total charge, atom id
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
else
    exit 1
fi
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
 
    #Generate copies for the files 
    cp COMPLEX.mol2 oldCOMPLEX.mol2
    cp COMPLEX.lib  oldCOMPLEX.lib

    # Reshape the charges in complex_esp.chg to be one column and save to charges.dat
    awk '{for (i=1; i<=NF; i++) print $i}' "$RUN_DIR/complex_esp.chg" > "$RUN_DIR/charges.chg"

    # Replace the charges in COMPLEX_charged.mol2 with those from charges.dat
    python3 "$SCRIPT_DIR/Retrieve_RESP_Charges.py"
    
    cp "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/COMPLEX_charged.mol2"
 
    tleap -f input_library.tleap > leap.log

    #Function to correct lib file
    python3 "$SCRIPT_DIR/12_generate_lib.py"
    
    #Generate the correct fixed and not fixed charges
    if [[ "${charmm_FF,,}" =~ ^(y|yes)$ ]] && [[ "${metalloprotein_choice,,}" =~ ^(n|no)$ ]] ; then
	    python3 "$SCRIPT_DIR/generate_str.py" "$RUN_DIR/COMPLEX.mol2" "$RUN_DIR/forcefield.dat" > "$RUN_DIR/COMPLEX_charged.top" 
	    echo " "
	    echo "Charges updated in COMPLEX_charged.top"
	    echo " "
    else
	    cp "$RUN_DIR/COMPLEX.lib" "$RUN_DIR/COMPLEX_charged.lib"
	    cp "$RUN_DIR/oldCOMPLEX.mol2" "$RUN_DIR/COMPLEX.mol2"
	    cp "$RUN_DIR/oldCOMPLEX.lib" "$RUN_DIR/COMPLEX.lib"
	    # Print output name
	    echo " "
	    echo "Charges updated in COMPLEX_charged.mol2"
	    echo "lib file updated in COMPLEX_charged.lib"
	    echo " "
    fi
	   

    # Clean up temporary files
    rm -f "$RUN_DIR/a" "$RUN_DIR/b" "$RUN_DIR/c" "$RUN_DIR/tmpfile" complex_esp.in complex_esp.out complex_esp.qin complex_esp.chg complex_esp.pch esp.dat charges.dat restraint_info.dat oldCOMPLEX.mol2 oldCOMPLEX.lib leap.log
else
    echo "No charge restraints applied. Exiting."
fi

files_to_remove=( "similar.dat" "input_library.tleap" "distance_type.dat" "tempz.fchk" "atomic_number.dat" "charges.dat" "resp.in" "temp.dat" "resp.out" "esp.chg" "bond_angle_dihedral_data.dat" "forcefield.dat" )

for file in "${files_to_remove[@]}"; do
    if [ -e "$file" ]; then
        rm -f "$file"
    fi
done
