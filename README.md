# 🚀 easyPARM v3.00 is Here!  

A new version of **easyPARM (v3.00)** has been released! 🎉  
This version includes **new features, enhancements, and expanded compatibility.**  

👉 **Check the latest release notes and manual **  


# **easyPARM v1.00**

**easyPARM** is a computational tool developed by **Abdelazim M. A. Abdelagawwad** under the supervision
of **Dr. Antonio Francés-Monerris** at the Institut de Ciència Molecular (ICMol), Universitat de València. 
Its primary purpose is to simplify the derivation of force field parameters for metal-containing molecular 
systems and to enable charge restraints on specific atoms.

## How to cite 
If you use this code in your research, please cite the original article associated with the code.

*Abdelgawwad AMA, Francés-Monerris A. easyPARM: Automated, Versatile, and Reliable Force Field Parameters for Metal-Containing Molecules with Unique Labeling of Coordinating Atoms. ChemRxiv. 2024;* https://doi.org/10.26434/chemrxiv-2024-f8wp4

## Key Features
1. **Force Field Parameter Generation**: easyPARM uses the Seminario method [1] to derive bond-stretching 
and bond-angle bending parameters for metal-containing systems. These parameters are optimized for use 
with the AMBER software suite [2], based on the Hessian matrix from frequency calculations performed with 
Gaussian 09 or 16. Additionally, the tool can generate parameters for non-metal systems using GAFF or the 
AMBER force fields. 

2. **Charge Restraint Using REsP Fitting**: easyPARM also supports charge restraints using the Restrained 
Electrostatic Potential (REsP) fitting approach [3]. This ensures that the electrostatic potential around 
the molecule is accurate and physically meaningful.

## Prerequisites
Before running easyPARM, you will need the following input files:

1. **Checkpoint (.chk) or Formatted Checkpoint File (.fchk)**: Contains the hessian matrix and data from a previous quantum 
chemistry calculation (Gaussian) [4].
2. **Optimized Structure (XYZ Format)**: A file with the optimized molecular geometry in XYZ format.
3. **Gaussian Charge Output**: Charge distribution information (ESP charge) generated from a Gaussian 
calculation.

## Running easyPARM
After downloading the code, if the files `easyPARM.sh`, `01_easyPARM.sh`, `04_parmch2_frcmod.sh`, and `readit` do not have executable permissions, please change their permissions using the following command:

`chmod +x file_name `

To run the code, you have two methods

**Method 1: Direct Script Execution**
To run the script directly, use the following command:

`./easyPARM.sh`

**Method 2: Using an Alias**
For easier access, you can create an alias in your .bashrc file:

`alias easyPARM='/full/path/to/easyPARM.sh'`

## Installation Requirements

**Python 3**: Ensure Python 3 is installed in your environment.

**Required Python Packages**: The script depends on scipy and periodictable. Install them using pip if not 
already available:

`pip install scipy` 

`pip install periodictable`

Additional package requirements may arise during execution. If prompted, install them via pip or contact 
us for assistance: abdelazim.abdelgawwad@uv.es.

## Important Notes
- Ensure that input files are correctly formatted and placed in the appropriate directory before running 
the script.
- The script processes these inputs to generate the necessary AMBER force field parameters, which will be 
output in the same directory.
- Output files will include `.mol2`, `.lib`, `.pdb`, and `.frcmod` files for use in molecular dynamics 
simulations.

## Molecular Dynamics Considerations
When using the mol2 and frcmod files to generate a library, tleap in AMBER may not correctly assign 
connectivity and atomic numbers for metals. To avoid this issue, we recommend using the pre-generated 
`.lib` and `.frcmod files`, which contain all necessary data.

Your leap script should include:


`loadamberparams COMPLEX.frcmod`

`loadoff COMPLEX.lib`

## Contact & Support
For troubleshooting or additional assistance, please refer to the documentation in our GitHub repository. 
You can also reach us via email at: abdelazim.abdelgawwad@uv.es.

## How to cite 
If you use this code in your research, please cite the original article associated with the code.


### References
[1]	Seminario, J. M. Calculation of Intramolecular Force Fields from Second-Derivative Tensors. Int J Quantum Chem 1996, 60 (7), 1271–1277. https://doi.org/https://doi.org/10.1002/(SICI)1097-461X(1996)60:7<1271::AID-QUA8>3.0.CO;2-W.

[2]	Case, D. A.; Aktulga, H. M.; Belfon, K.; Ben-Shalom, I. Y.; Brozell, S. R.; Cerutti, D. S.; Cheatham 
III, T. E.; Cisneros, G. A.; Cruzeiro, V. W. D.; Darden, T. A.; Duke, R. E.; Giambasu, G.; Gilson, M. K.; 
Gohlke, H.; Goetz, A. W.; Harris, R.; Izadi, S.; Izmailov, S. A.; Jin, C.; Kasavajhala, K.; Kaymak, M. C.; King, 
E.; Kovalenko, A.; Kurtzman, T.; Lee, T. S.; LeGrand, S.; Li, P.; Lin, C.; Liu, J.; Luchko, T.; Luo, R.; 
Machado, M.; Man, V.; Manathunga, M.; Merz, K. M.; Miao, Y.; Mikhailovskii, O.; Monard, G.; Nguyen, H.; 
Oâ€™Hearn, K. A.; Onufriev, A.; Pan, F.; Pantano, S.; Qi, R.; Rahnamoun, A.; D.R. Roe; Roitberg, A.; Sagui, C.; 
Schott-Verdugo, S.; Shen, J.; Simmerling, C. L.; Skrynnikov, N. R.; Smith, J.; Swails, J.; Walker, R. C.; Wang, 
J.; Wei, H.; Wolf, R. M.; Wu, X.; Xue, Y.; York, D. M.; Zhao, S.; Kollman, P. A. Amber 20. Amber 20 2021, San 
Francisco.      https://ambermd.org/

[3]	Comell, W. D.; Cieplak, P.; Bayly, C. I.; Kollman, P. A. Application of RESP Charges to calculate 
conformational energies, hydrogen bond energies, and free energies of solvation. Journal of American Chemical 
Society 1993, 115 (7), 9620-9631.   https://doi.org/10.1021/ja00074a030

[4] Frisch, M. J.; Trucks, G. W.; Schlegel, H. B.; Scuseria, G. E.; Robb, M. a.; Cheeseman, J. R.; Scalmani, 
G.; Barone, V.; Petersson, G. a.; Nakatsuji, H.; Li, X.; Caricato, M.; Marenich, a. V.; Bloino, J.; Janesko, 
B. G.; Gomperts, R.; Mennucci, B.; Hratchian, H. P.; Ortiz, J. V.; Izmaylov, a. F.; Sonnenberg, J. L.; 
Williams; Ding, F.; Lipparini, F.; Egidi, F.; Goings, J.; Peng, B.; Petrone, A.; Henderson, T.; Ranasinghe, 
D.; Zakrzewski, V. G.; Gao, J.; Rega, N.; Zheng, G.; Liang, W.; Hada, M.; Ehara, M.; Toyota, K.; Fukuda, R.; 
Hasegawa, J.; Ishida, M.; Nakajima, T.; Honda, Y.; Kitao, O.; Nakai, H.; Vreven, T.; Throssell, K.; Montgomery 
Jr., J. a.; Peralta, J. E.; Ogliaro, F.; Bearpark, M. J.; Heyd, J. J.; Brothers, E. N.; Kudin, K. N.; 
Staroverov, V. N.; Keith, T. a.; Kobayashi, R.; Normand, J.; Raghavachari, K.; Rendell, a. P.; Burant, J. C.; 
Iyengar, S. S.; Tomasi, J.; Cossi, M.; Millam, J. M.; Klene, M.; Adamo, C.; Cammi, R.; Ochterski, J. W.; 
Martin, R. L.; Morokuma, K.; Farkas, O.; Foresman, J. B.; Fox, D. J. G16_C01. 2016, p Gaussian 16, Revision C.01, 
Gaussian, Inc., Wallinford, CT, USA.       https://gaussian.com/ 

# easyPARM v2.00 Release Notes
## What's New

- **ORCA Output Parsing**: Version 2.00 adds support for parsing ORCA [1] output files, enabling users to process ORCA data along with previously supported Gaussian outputs.
- **Atomic Charges with CHELPG**: Now supports atomic charge calculations using the CHELPG scheme [2] from ORCA.
- **Enhanced Complex Compatibility**: Works with a wider range of non-standard complexes, including those with inorganic ligands like COSAN.

## Running Scheme Update

When running easyPARM, you now have the following options to choose from, depending on the charge output and QM output type you are working with:

### Charge Calculation Method
`Select the charge calculation method:`

`1- Gaussian (RESP charges)`

`2- ORCA (CHELPG charges)`

### QM Output Format
`Please select the format you will provide:`

`1- Orca Output`

`2- Gaussian Output`

`3- Gaussian Checkpoint`

`4- Gaussian Formatted Checkpoint`

After selecting your desired option, provide the corresponding required files as described below.

## Required files
- **Orca Output:** For Orca, you must provide both the standard output file and the Hessian output (output.hess). These files are generated by running a **single-point frequency calculation.**

- **Gaussian Output:** For Gaussian, you only need the output file. However, ensure that you run a single-point frequency calculation and include the keyword `iop(7/33=1)` in your input file.

- **Checkpoint and Formatted Checkpoint:** These options remain unchanged from the previous version.

## Additional Requirements

Regardless of the selected output type, make sure to also provide:

- **The charge output**
- **The optimized geometry in XYZ format**

## Molecular Dynamics Considerations
### Important Note

If your system contains a metal with a coordination bond greater than 8, you may encounter the following warning when generating the **prmtop** file in `tleap`:

`Bond: Maximum coordination exceeded on A<Co1 23>`
     
   `-- setting atoms pert=true overrides default limits`
      
This warning indicates that the metal has more bonds than the allowed limit. Although tleap will generate the prmtop file, it will be incorrect because not all bonds for the metal will be included.

### Solution
To address this, modify the allowed bond limit in the AmberTools source code:
1. Access the AmberTools path:
`$AMBER_HOME/AmberTools/src/leap/src/leap`
2. Edit the atom.h file to increase the bond limit. By default, the maximum number of bonds is set to 8:

`/* maximum 8 bonds out of each atom */`

`#define MAXBONDS 8`

3. Change 8 to the desired number of bonds.
4. Recompile Amber to apply the changes, allowing tleap to recognize the new MAXBONDS limit.

## References
1- Frank Neese, Frank Wennmohs, Ute Becker, Christoph Riplinger; The ORCA quantum chemistry program package. J. Chem. Phys. 14 June 2020; 152 (22): 224108. https://doi.org/10.1063/5.0004608
2- Breneman, C.M. and Wiberg, K.B. (1990), Determining atom-centered monopoles from molecular electrostatic potentials. The need for high sampling density in formamide conformational analysis. J. Comput. Chem., 11: 361-373. https://doi.org/10.1002/jcc.540110311

# easyPARM v3.00 Release Notes  

## 🚀 What's New in v3.00  

- **Force Field Parameterization for Metalloproteins**  
- **Expanded Charge Calculation Support**  
  - RESP charges via **ORCA**  
  - RESP and other charge models via **GAMESS**  
- **Amber Format Conversion**  
  - Export to **GROMACS** or **OpenMM**  

📖 **Check the manual for more details.**  

Note:- For metalloprotein, please, *install Biopython library*.


