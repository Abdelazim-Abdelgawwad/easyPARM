#!/bin/bash

# Capture the directory from which this script is called
ORIGINAL_DIR="$(pwd)"

# Determine the directory of this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Navigate to the scripts directory (if needed)
cd "$SCRIPT_DIR/scripts"

# Run the 01_easyPARM.sh script, passing the original working directory as an argument
./01_easyPARM.sh "$ORIGINAL_DIR"

