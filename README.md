# 🚀 easyPARM v4.15 Released! 🎉

We’re excited to announce the release of **easyPARM v4.15**, featuring several **bug fixes**, **user-focused debugging improvements**, and an **expanded tutorial section**!

## 🛠️ Bug Fixes & Improvements
- Fixed issues related to **Gaussian** and **ORCA** output parsing.  
- Resolved problems affecting specific **atom types**.  
- Improved **debug messages** for clearer error reporting to users.  
- Added a **new tutorial** dedicated to **metalloprotein parameterization**.

## 🔥 Previoulsy in v4.10

🧪 **Psi4 Integration**
easyPARM now supports the **[Psi4](https://psicode.org/psi4manual/master/index.html)** open-source quantum chemistry package.
All you need is your input `.xyz` structure — easyPARM v4.10 will handle:

* **Geometry optimization**
* **Frequency analysis**
* **Charge calculations**
* **Force field parameter generation**

📖 Check the detailed [Psi4 Support Documentation](https://abdelazim-abdelgawwad.github.io/Tutorial/Psi4%20Support).

🔗 Within the Psi4 framework, easyPARM uses the [RESP implementation by Alenaizan](https://github.com/cdsgroup/resp).

⚡ **Enhanced ORCA Support**  
- Older versions of **ORCA** did not include RESP fitting.  
- With **ORCA 6.0**, the new `!RESP` keyword enables RESP charge fitting directly within ORCA’s scheme.  
- In addition, **easyPARM** provides its own RESP fitting procedure using user-specified electrostatic potential (ESP) points.  
- This means you can now choose between **ORCA’s built-in RESP** or **easyPARM’s RESP workflow**, depending on your needs.  

🔧 **Bug Fixes & Stability**

* Numerous user-reported bugs fixed — special thanks to the community for feedback and contributions!
* Enhanced handling of **symmetry** and **equivalent atom environments** for more robust parameter generation.

---

## 🔄 Previously in v4.00

🧠 **Improved Atom Type Detection**

* Smarter recognition of **challenging coordination complexes**.
* Accurate typing via **molecular topology, hybridization, ring strain, aromaticity, and metal–ligand recognition**.

🧬 **ff19SB for Metal-Coordinated Residues**

* Standard residues bound to metals now typed using **ff19SB**, improving consistency in **metalloproteins**.

🔬 **Updated Non-Bonded Parameters for Metals**

* Monovalent ions from [Li et al., *JCIM* (2021)](https://doi.org/10.1021/acs.jcim.0c01390)
* Divalent ions from [Zhang et al., *JCTC* (2021)](https://doi.org/10.1021/acs.jctc.0c00194)
* Ensures full **AMBER compatibility**.

🚨 **Bug Fixes** for improved stability.

---

## 📚 Documentation

📖 [Tutorial & Documentation](https://abdelazim-abdelgawwad.github.io/Tutorial/)
📂 Explore usage examples, Psi4 workflows, and advanced configuration tips.

---

## 💡 Contribute & Suggest Improvements

We welcome **issues**, **pull requests**, and **discussions**:

* [Open an Issue](https://github.com/Abdelazim-Abdelgawwad/easyPARM/issues)
* [Submit a Pull Request](https://github.com/Abdelazim-Abdelgawwad/easyPARM/pulls)
* Join the community in the **Discussions** tab

Let’s build the future of **automated molecular parameterization** together! 🧪⚡

---

## 📖 How to Cite

If you use **easyPARM** in your research, please cite:

> Abdelazim M. A. Abdelgawwad and Antonio Francés-Monerris.
> *easyPARM: Automated, Versatile, and Reliable Force Field Parameters for Metal-Containing Molecules with Unique Labeling of Coordinating Atoms.*
> *J. Chem. Theory Comput. 2025, 21, 4, 1817–1830*
> [DOI: 10.1021/acs.jctc.4c01272](https://doi.org/10.1021/acs.jctc.4c01272)

---
