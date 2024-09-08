'easyPARM'

Instructions for Running easyPARM
To run the easyPARM code, you will need to provide the following three input files:

1-Checkpoint File (.chk): This file contains the wavefunction and other necessary data from a previous quantum chemistry calculation (Gaussian).
2-Optimized Structure (XYZ Format): A file containing the optimized molecular geometry in XYZ format.
3-Gaussian Charge Output: A file generated from a Gaussian calculation, containing the charge distribution information (esp charge).

Running the Code
Once you have the required files, you can run the code using the following command:
./easyPARM.sh

Important Notes
Python 3: Please ensure that Python 3 is installed and available in your environment. The code requires Python 3 to function correctly.

Required Python Package: The code depends on the periodictable Python package. If this package is not already installed, you can install it using pip: pip install periodictable

Additional Information
Before running the code, make sure your input files are correctly formatted and placed in the appropriate directory.
The script is designed to automatically process the provided input files and generate the necessary AMBER force field parameters.
Output files will be generated in the same directory where the script is executed, containing the newly created force field parameters ready for use in molecular dynamics simulations.

For troubleshooting or further assistance, please refer to the documentation available on our GitHub repository.
By following these steps, you can easily generate accurate force field parameters for metal-containing systems using easyPARM.
