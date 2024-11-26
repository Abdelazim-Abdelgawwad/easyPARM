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
#                              |  $$$$$$/              Ver. 2.00 - 1 November 2024                                #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import sys
from collections import defaultdict

def xyz_to_pdb(xyz_file, pdb_file):
    with open(xyz_file, 'r') as xyz, open(pdb_file, 'w') as pdb:
        # Skip the first two lines of the XYZ file
        next(xyz)
        next(xyz)
        
        atom_number = 1
        atom_counts = defaultdict(int)
        
        for line in xyz:
            atom, x, y, z = line.split()
            atom_counts[atom] += 1
            atom_label = f"{atom}{atom_counts[atom]}"
            
            pdb.write(f"ATOM  {atom_number:5d} {atom_label:<3} mol      1    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}  1.00  0.00          \n")
            atom_number += 1
        
        pdb.write("END\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    pdb_file = sys.argv[2]
    
    xyz_to_pdb(xyz_file, pdb_file)
