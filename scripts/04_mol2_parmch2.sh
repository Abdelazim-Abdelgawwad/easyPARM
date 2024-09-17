#!/bin/bash

# 1. Ask the user for the path to Amber

# Get the directory of the bash script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Remove the existing COMPLEX_modified.frcmod file if it exists
file_name="COMPLEX_modified.frcmod"

if [ -f "$file_name" ]; then
    rm "$file_name"
fi

# Run the Python script using the directory of the bash script
python3 "$SCRIPT_DIR/generate_preforcefield.py"

# Run parmchk2 with the appropriate parameters
parmchk2 -i COMPLEX_modified.mol2 -f mol2 -o COMPLEX_modified.frcmod -a Y

# Run the Python script again
python3 "$SCRIPT_DIR/generate_preforcefield.py"

