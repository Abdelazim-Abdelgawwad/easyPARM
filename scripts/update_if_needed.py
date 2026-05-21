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
#                              |  $$$$$$/              Ver. 4.20 - 1 January 2026                                 #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

import re
import os
import sys
import shutil
BOND_SCALE_RULES = [
    (20, 4.599),   
]
ANGLE_SCALE_RULES = [
    (5,  11.599),  
    (10,  7.799),  
    (20,  3.599),  
    (29,  2.699),  
]
LINE_RE = re.compile(
    r"^([A-Za-z0-9*']+"
    r"(?:\s*-\s*[A-Za-z0-9*']+)+"  # one or more  -ATOM segments
    r")"
    r"(\s+)([-0-9.]+)"              # group 3 = force
    r"(\s+)([-0-9.]+)"             # group 5 = eq_val
    r"(.*)"                         # group 6 = trailing comment (preserved)
)
#Return scaled force value according to the first matching rule, or unchanged.
def apply_scale(force, rules):
    for threshold, factor in rules:
        if force < threshold:
            return force * factor
    return force  # no rule matched → unchanged
def process_frcmod(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    section = None
    output = []
    updated = 0
    unchanged = 0
    for line in lines:
        stripped = line.rstrip('\n')
        if re.match(r'^BOND', stripped):
            section = 'BOND'
            output.append(line)
            continue
        elif re.match(r'^ANGLE', stripped):
            section = 'ANGLE'
            output.append(line)
            continue
        elif re.match(r'^DIHE|^IMPROPER|^NONBON', stripped):
            section = stripped[:6]
            output.append(line)
            continue
        if not stripped or stripped.startswith('#'):
            output.append(line)
            continue
        if section in ('BOND', 'ANGLE') and re.match(r'^[A-Za-z]', stripped):
            m = LINE_RE.match(stripped)
            if m:
                atom_pair  = m.group(1)   # e.g. "C -O3" or "CA-N"
                sep1       = m.group(2)   # whitespace before force
                force_str  = m.group(3)   # e.g. "622.90"
                sep2       = m.group(4)   # whitespace before eq_val
                eq_val     = m.group(5)   # e.g. "1.225"
                trailing   = m.group(6)   # any trailing comment
                force = float(force_str)
                if section == 'BOND':
                    new_force = apply_scale(force, BOND_SCALE_RULES)
                    fmt = f"{new_force:.2f}"
                else:
                    new_force = apply_scale(force, ANGLE_SCALE_RULES)
                    fmt = f"{new_force:.2f}"
                new_line = f"{atom_pair}{sep1}{fmt}{sep2}{eq_val}{trailing}\n"
                output.append(new_line)
                if new_force != force:
                    msg = "  [%s] %-20s  force: %s -> %s" % (section, atom_pair.strip(), force_str, fmt)
                    updated += 1
                else:
                    unchanged += 1
                continue
        output.append(line)
    backup = filepath + ".bak"
    shutil.copy2(filepath, backup)
    with open(filepath, 'w') as f:
        f.writelines(output)
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 scale_if_needed.py <path/to/COMPLEX.frcmod>")
        sys.exit(1)
    frcmod_path = sys.argv[1]
    if not os.path.exists(frcmod_path):
        print(f"Error: File not found: {frcmod_path}")
        sys.exit(1)
    process_frcmod(frcmod_path)
if __name__ == "__main__":
    main()

