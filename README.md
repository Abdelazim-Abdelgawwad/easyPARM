# **easyPARM**

**easyPARM** is a computational tool developed by **Abdelazim M. A. Abdelagawwad** under the supervision
of **Dr. Antonio Francés-Monerris** at the Institut de Ciència Molecular (ICMol), Universitat de València. 
Its primary purpose is to simplify the derivation of force field parameters for metal-containing molecular 
systems and to enable charge restraints on specific atoms.

## How to cite 
If you use this code in your research, please cite the original article associated with the code.

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
[1]	Seminar, J. M. Calculation of Intramolecular Force Fields from Second-Derivative Tensors; 1996. https://doi.org/https://doi.org/10.1002/(SICI)1097-461X(1996)60:7<1271::AID-QUA8>3.0.CO;2-W

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

