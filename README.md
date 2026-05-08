# easyPARM

![Views](https://hits.sh/github.com/Abdelazim-Abdelgawwad/easyPARM.svg?style=for-the-badge&label=Views&color=2563eb&labelColor=1e3a8a)
[![GitHub stars](https://img.shields.io/github/stars/Abdelazim-Abdelgawwad/easyPARM?style=for-the-badge&logo=github&color=f59e0b&labelColor=78350f)](https://github.com/Abdelazim-Abdelgawwad/easyPARM/stargazers)
![GitHub forks](https://img.shields.io/github/forks/Abdelazim-Abdelgawwad/easyPARM?style=for-the-badge&color=f97316&labelColor=7c2d12)

[![Citations JCTC 2025](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.4c01272&query=%24.message&style=for-the-badge&label=Citations%20JCTC%202024&color=16a34a&labelColor=14532d&logo=google-scholar&logoColor=white)](https://doi.org/10.1021/acs.jctc.4c01272)

[![Citations JCP 2025](https://img.shields.io/badge/Citations%20JCP%202025-0-7c3aed?style=for-the-badge&logo=google-scholar&logoColor=white&labelColor=4c1d95)](https://doi.org/10.1063/5.0301038)

----

# 🚀 New Release: easyPARM v4.20

We are pleased to announce **easyPARM v4.20**, a major update that significantly expands support for **metal–nucleic acid systems**, enhances **metalloprotein force-field flexibility**, and resolves several critical bugs in the parameterization workflow.

---

## 🧬 Metal–Nucleic Acid Support (New in v4.20)

easyPARM now **fully supports metal–nucleic acid systems**, enabling robust parameterization of a wide range of biologically and technologically relevant constructs, including:

* **Metal-bound DNA and RNA**
* **Inner-sphere metal coordination** to nucleobases and phosphate groups
* **Metal-mediated base pairs** and metallo-DNA architectures
* **Nucleic acid–metal complexes** and nucleic acid–directed metal assemblies

### 🔧 Force Field Selection for Coordinated Nucleotides

Users can now explicitly select the force field for **metal-coordinated standard nucleotides**, including:

* **OL15**
* **OL21**
* **OL24**
* **BSC0**
* **BSC1**
* **OL3**

This provides full flexibility when working with modern AMBER nucleic acid force fields in metal-coordination environments.

📖 A dedicated tutorial with practical examples is available here:
[https://abdelazim-abdelgawwad.github.io/Tutorial/Tutorial/metallonucleicacid](https://abdelazim-abdelgawwad.github.io/Tutorial/Tutorial/metallonucleicacid)

---

## 🧪 Enhanced Metalloprotein Force-Field Support

Previously, metalloprotein parameterization in easyPARM required the exclusive use of **ff19SB**.
With **v4.20**, users can now choose among multiple AMBER protein force fields:

* **ff12SB**
* **ff14SB**
* **ff19SB**
* **AMBER-fb15**

This update enables greater methodological consistency across mixed systems and legacy workflows.

---

## 🛠️ Bug Fixes & Stability Improvements

### Metalloprotein Parameterization

Several issues affecting metalloprotein workflows have been resolved, including:

1. **Incorrect bond information** in the `Bond_info.dat` file
2. **Improper handling of terminal residues** when coordinated to metal centers

### General Fixes

* Fixed bugs affecting the **generation of certain MOL2 files**
* Improved robustness and consistency of intermediate parameter files

---

## 🔄 Previously: easyPARM v4.15

**easyPARM v4.15** introduced important stability improvements and usability enhancements:

### 🛠️ Bug Fixes & Debugging

* Fixed issues in **Gaussian** and **ORCA** output parsing
* Resolved problems with specific **atom-type assignments**
* Improved **debug messages** for clearer user feedback

### 📖 Documentation

* Added a **dedicated metalloprotein parameterization tutorial**

---

## 🔥 Previously: easyPARM v4.10

### 🧪 Psi4 Integration

Support for the **Psi4** open-source quantum chemistry package, enabling automated:

* Geometry optimization
* Frequency calculations
* RESP charge derivation
* Force-field parameter generation

📖 Psi4 documentation:
[https://abdelazim-abdelgawwad.github.io/Tutorial/Psi4%20Support](https://abdelazim-abdelgawwad.github.io/Tutorial/Psi4%20Support)

easyPARM integrates the **RESP implementation by Alenaizan** within the Psi4 framework.

### ⚡ Enhanced ORCA Support

* Native RESP fitting via `!RESP` in **ORCA 6.0**
* Optional **easyPARM RESP workflow** using user-defined ESP points
* Users can freely choose between **ORCA-native** or **easyPARM-based** RESP fitting

---

## 🧠 Previously: easyPARM v4.00

* Improved **atom-type detection** for challenging coordination environments
* Accurate recognition based on **topology, hybridization, aromaticity, ring strain, and metal–ligand bonding**
* Introduction of **ff19SB typing** for metal-coordinated standard residues
* Updated **non-bonded metal parameters** compatible with AMBER:

  * Monovalent ions: Li *et al.*, *JCIM* (2021)
  * Divalent ions: Zhang *et al.*, *JCTC* (2021)

---

## 📚 Documentation & Tutorials

📖 Full documentation and tutorials:
[https://abdelazim-abdelgawwad.github.io/Tutorial/](https://abdelazim-abdelgawwad.github.io/Tutorial/)

---

## 💡 Community & Contributions

Feedback and contributions are welcome:

* 🐛 Issues: [https://github.com/Abdelazim-Abdelgawwad/easyPARM/issues](https://github.com/Abdelazim-Abdelgawwad/easyPARM/issues)
* 🔧 Pull Requests: [https://github.com/Abdelazim-Abdelgawwad/easyPARM/pulls](https://github.com/Abdelazim-Abdelgawwad/easyPARM/pulls)
* 💬 Discussions: GitHub Discussions tab

---

## 📖 How to Cite

If you use **easyPARM** in your research, please cite:

> Abdelazim M. A. Abdelgawwad and Antonio Francés-Monerris
> *easyPARM: Automated, Versatile, and Reliable Force Field Parameters for Metal-Containing Molecules with Unique Labeling of Coordinating Atoms*
> **J. Chem. Theory Comput.** 2025, 21, 4, 1817–1830
> [DOI: 10.1021/acs.jctc.4c01272](https://doi.org/10.1021/acs.jctc.4c01272)

> Abdelazim M. A. Abdelgawwad and Antonio Francés-Monerris.
> *easyPARM v4.00: A Python-based tool for the automated parameterization of metalloproteins and metal–organic polyhedra with multiple metal centers.*
> *Journal of Chemical Physics*, 2025, **163**(22), 222501.
> [DOI: 10.1063/5.0301038](https://doi.org/10.1063/5.0301038)
