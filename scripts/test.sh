#!/bin/bash

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

echo " "
echo " "
read -p "Would you like to restrain the charge on specific atoms? (Yes or No): " restrain_choice
if [[ "$restrain_choice" =~ ^[Yy][Ee][Ss]|[Yy]$ ]]; then
    echo " "
    while true; do
        num_atoms=$(get_valid_input "How many atoms do you want to restrain? " "1 2")
        
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

        echo "Restraint information has been saved to $RUN_DIR/restraint_info.dat"
        
        # Ask if the user wants to review or redo the restraints
        read -p "Do you want to review the restraints? (y/n): " review_choice
        if [[ "$review_choice" =~ ^[Yy]$ ]]; then
            cat "$RUN_DIR/restraint_info.dat"
            read -p "Do you want to redo the restraints? (y/n): " redo_choice
            if [[ "$redo_choice" =~ ^[Yy]$ ]]; then
                continue
            fi
        fi

        break
    done
    
    echo " "
    read -p "Please enter the name of Gaussian output that contains all the ESP charges: " ESP_FILE

    # Extract relevant data from the Gaussian output
    grep "Atomic Center " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/a"
    grep "ESP Fit" "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/b"
    grep "Fit    " "$RUN_DIR/$ESP_FILE" > "$RUN_DIR/c"

    # Execute the binary that processes this data
    i=$(wc -l < "$RUN_DIR/a")
    
    # Ensure the readit binary is in the script directory and is executable
    if [ ! -x "$SCRIPT_DIR/readit" ]; then
        echo "Error: readit binary not found or not executable in $SCRIPT_DIR"
        exit 1
    fi
    
    # Run the readit binary
    (cd "$RUN_DIR" && "$SCRIPT_DIR/readit")
    

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

# Read the atom numbers and charges to restrain into an associative array
    declare -A restrained_atoms
    while read -r atom_number charge; do
    	restrained_atoms[$atom_number]=$charge
    done < "$RUN_DIR/restraint_info.dat"

# Read the similar.dat file into another associative array for atom types
    declare -A atom_types
    while read -r atom_number atom_type; do
    	atom_types[$atom_number]=$atom_type
    done < "$RUN_DIR/similar.dat"

# Convert associative arrays to space-separated lists
    restrained_atoms_list=$(for key in "${!restrained_atoms[@]}"; do echo "$key ${restrained_atoms[$key]}"; done)
    atom_types_list=$(for key in "${!atom_types[@]}"; do echo "$key ${atom_types[$key]}"; done)

# Extract and process the data
    awk -v i="$i" -v OFILE="$OUTPUT_FILE" -v restrained_atoms_list="$restrained_atoms_list" -v atom_types_list="$atom_types_list" '
BEGIN {
    start=0;
    count=0;
    
    # Populate associative arrays in awk
    split(restrained_atoms_list, restrained_array, " ");
    for (i = 1; i <= length(restrained_array); i+=2) {
        restrained_atoms[restrained_array[i]] = restrained_array[i+1];
    }
    
    split(atom_types_list, atom_types_array, " ");
    for (i = 1; i <= length(atom_types_array); i+=2) {
        atom_types[atom_types_array[i]] = atom_types_array[i+1];
    }
}
    
/Standard orientation:/ {start=1; next}

start && /^ *[0-9]+ +[0-9]+ +[0-9]+ +[-.0-9]+ +[-.0-9]+ +[-.0-9]+$/ {
    count++;
    if (count <= i) {
        atom_num = $1;
        atom_type = $3;
        
        if (atom_num in restrained_atoms) {
            atom_type = -1;
        } else if (atom_num in atom_types) {
            atom_type = atom_types[atom_num];
        }
        
        printf "%3d%4d\n", $2, atom_type >> OFILE;
    }
}

/^---------------------------------------------------------------------$/ {start=0; count=0}
' "$RUN_DIR/$ESP_FILE"
 
    # File to read
    echo " " >> complex_esp.in
    echo " " >> complex_esp.in

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

    echo " "
    echo "Charges updated in COMPLEX_charged.mol2."

    # Clean up temporary files
    rm -f "$RUN_DIR/a" "$RUN_DIR/b" "$RUN_DIR/c" "$RUN_DIR/tmpfile" complex_esp.in complex_esp.out complex_esp.qin complex_esp.chg complex_esp.pch esp.dat charges.dat restraint_info.dat  
else
    echo "No charge restraints applied. Exiting."
fi

rm similar.dat

