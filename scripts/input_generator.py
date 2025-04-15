#!/usr/bin/env python3
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
#                              |  $$$$$$/              Ver. 3.25 - 14 April 2025                                  #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################


import sys
import re

#Check the input file and return a dictionary of parameters.
def check_input_file(input_path):
    parameters = {}
    
    try:
        with open(input_path, 'r') as f:
            for line in f:
                # Remove comments (anything after #)
                line = line.split('#', 1)[0].strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Read key-value pairs (KEY = VALUE format)
                match = re.match(r'^([A-Za-z0-9_$]+)\s*=\s*(.+)$', line)
                if match:
                    key = match.group(1).upper()
                    value = match.group(2).strip()
                    parameters[key] = value
    except Exception as e:
        print(f"Error parsing input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return parameters

#Generate the input sequence for 01_easyPARM.sh based on parameters.
def generate_input_for_easyparm(parameters):
    inputs = []
    missing_params = []
    
    # Helper function to get parameter or track if missing
    def get_param(key, default=None, required=True):
        if key in parameters:
            return parameters[key]
        elif required:
            missing_params.append(key)
            return default
        else:
            return default
    
    # 1. PARAMETERIZE
    parameterize_value = get_param("PARAMETERIZE", "1", required=False)
    
    # 2. AMBER_CONFIG
    amber_config = get_param("AMBER_CONFIG", "1")
    inputs.append(amber_config)
    
    # Add AMBER_PATH if AMBER_CONFIG is 2
    if amber_config == "2":
        inputs.append(get_param("AMBER_PATH", "/usr/local/amber"))

    # 3. CHARGE
    inputs.append(get_param("CHARGE", "0"))
    
    # 4. MULTIPLICITY
    inputs.append(get_param("MULTIPLICITY", "1"))
    
    # 5. STRUCTURE_FILE
    inputs.append(get_param("STRUCTURE_FILE", "OPTIMIZED.xyz"))
    
    # 6. CHARGE_METHOD
    charge_method = get_param("CHARGE_METHOD", "1")
    inputs.append(charge_method)
    
    # Conditional inputs based on CHARGE_METHOD
    if charge_method == "1":  # RESP (GAUSSIAN)
        inputs.append(get_param("INPUT_FORMAT", "1"))
        inputs.append(get_param("METHOD", "1"))
        inputs.append(get_param("ATOM_TYPE", "1"))
        inputs.append(get_param("CHARGE_OUTPUT", "output.log"))
    elif charge_method == "2":  # CHELPG (ORCA)
        inputs.append(get_param("ATOM_TYPE", "1"))
        inputs.append(get_param("CHARGE_OUTPUT", "output.out"))
    elif charge_method == "3":  # RESP (ORCA)
        inputs.append(get_param("ATOM_TYPE", "1"))
        inputs.append(get_param("CHARGE_OUTPUT", "output.vpot"))
    elif charge_method == "4":  # RESP (GAMESS)
        inputs.append(get_param("ATOM_TYPE", "1"))
        inputs.append(get_param("CHARGE_OUTPUT", "output.dat"))
    elif charge_method == "5":  # GAMESS CHARGES
        inputs.append(get_param("ATOM_TYPE", "1"))
        inputs.append(get_param("CHARGE_OUTPUT", "output.log"))
    
    # 7. QM_OUTPUT
    qm_output = get_param("QM_OUTPUT", "2")
    inputs.append(qm_output)
    
    # Conditional inputs based on QM_OUTPUT
    if qm_output == "1":  # Orca
        if charge_method in ["1", "3", "4", "5"]:
            inputs.append(get_param("ORCA_OUTPUT", "output.out"))
            inputs.append(get_param("ORCA_HESS", "output.hess"))
        elif charge_method in ["2"]:
            inputs.append(get_param("ORCA_HESS", "output.hess"))
    elif qm_output == "2":  # Gaussian Output
        inputs.append(get_param("GAUSSIAN_OUTPUT", "output.log"))
    elif qm_output == "3":  # Gaussian Checkpoint
        inputs.append(get_param("CHK_FILE", "output.chk"))
        gaussian_conf = get_param("GAUSSIAN_CONFIG", "1")
        inputs.append(gaussian_conf)
        if gaussian_conf == "2":
            inputs.append(get_param("GAUSSIAN_PATH", "/path/to/gaussian"))
    elif qm_output == "4":  # Gaussian Formatted Checkpoint
        inputs.append(get_param("FCHK_FILE", "output.fchk"))
    elif qm_output == "5":  # Gamess Output
        inputs.append(get_param("GAMESS_OUTPUT", "output.dat"))
    
    # 8. METALLOPROTEIN
    metalloprotein = get_param("METALLOPROTEIN", "no", required=False).lower()
    if metalloprotein in ['y', 'yes', 'true', '1']:
        inputs.append('y')
        # Include PDB file if metalloprotein is yes
        inputs.append(get_param("PDB_FILE", "protein.pdb"))
    else:
        inputs.append('n')
        
    # 9. RESID_ID
    resid_id = get_param("RESID_ID", "no", required=False).lower()
    if resid_id in ['y', 'yes', 'true', '1']:
        inputs.append('y')
        inputs.append(get_param("RESID_NAME", "LIG"))
    else:
        inputs.append('n')
    
    # 10. CHARMM_FF
    charmm_ff = get_param("CHARMM_FF", "no", required=False).lower()
    if charmm_ff in ['y', 'yes', 'true', '1']:
        inputs.append('y')
    else:
        inputs.append('n')
    
    # 11. CHARGE_RESTRAIN
    charge_restrain = get_param("CHARGE_RESTRAIN", "no", required=False).lower()
    if charge_restrain in ['y', 'yes', 'true', '1']:
        inputs.append('y')
        
        # Get number of atoms to restrain
        num_atoms_str = get_param("NUM_ATOMS", "0")
        try:
            num_atoms = int(num_atoms_str)
            inputs.append(str(num_atoms))
            
            # Process each atom's restraint information
            for i in range(1, num_atoms + 1):
                atom_key = f"ATOM{i}"
                atom_value = get_param(atom_key, f"1 0.0")
                inputs.append(atom_value)
        except ValueError:
            missing_params.append("NUM_ATOMS (invalid integer value)")
            inputs.append("0")  # Default if num_atoms is invalid
            
        # 12. REVIEW_CHARGES (always 'n')
        inputs.append('n')
    
        # 13. CHARGE_FILE
        inputs.append(get_param("CHARGE_FILE", "output.vpot"))
    else:
        inputs.append('n')
    
    # Return missing parameters and inputs
    return missing_params, inputs

def main():
    if len(sys.argv) != 3:
        print("Usage: input_generator.py input_file output_file", file=sys.stderr)
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    # Parse the input file
    parameters = check_input_file(input_path)
    
    # Generate input sequence for 01_easyPARM.sh
    missing_params, inputs = generate_input_for_easyparm(parameters)
    
    # Check if there are missing parameters
    if missing_params:
        print("\nError: Required parameters are missing from the input file:", file=sys.stderr)
        for param in missing_params:
            print(f"  - Missing parameter: '{param}'", file=sys.stderr)
        print("\nPlease correct your input file or continue with interactive mode", file=sys.stderr)
        sys.exit(1)
    
    # Write to output file only if no parameters are missing
    try:
        with open(output_path, 'w') as f:
            for input_value in inputs:
                f.write(f"{input_value}\n")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    sys.exit(0)

if __name__ == "__main__":
    main()
