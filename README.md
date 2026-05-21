# easyPARM

![Views](https://hits.sh/github.com/Abdelazim-Abdelgawwad/easyPARM.svg?style=for-the-badge\&label=Views\&color=2563eb\&labelColor=1e3a8a)
[![GitHub stars](https://img.shields.io/github/stars/Abdelazim-Abdelgawwad/easyPARM?style=for-the-badge\&logo=github\&color=f59e0b\&labelColor=78350f)](https://github.com/Abdelazim-Abdelgawwad/easyPARM/stargazers)
![GitHub forks](https://img.shields.io/github/forks/Abdelazim-Abdelgawwad/easyPARM?style=for-the-badge\&color=f97316\&labelColor=7c2d12)

[![Citations JCTC 2025](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.4c01272\&query=%24.message\&style=for-the-badge\&label=Citations%20JCTC%202024\&color=16a34a\&labelColor=14532d\&logo=google-scholar\&logoColor=white)](https://doi.org/10.1021/acs.jctc.4c01272)

[![Citations JCP 2025](https://img.shields.io/badge/Citations%20JCP%202025-0-7c3aed?style=for-the-badge\&logo=google-scholar\&logoColor=white\&labelColor=4c1d95)](https://doi.org/10.1063/5.0301038)

---

# 🚀 New Release: easyPARM v4.25

We are pleased to announce **easyPARM v4.25**, which introduces important improvements for **metalloprotein parameterization workflows**, enhanced automation for **tleap input generation**, and several bug fixes related to complex coordination environments.

---

## 🧪 Improved Metalloprotein Support

easyPARM v4.25 now correctly handles challenging metalloprotein coordination patterns where:

* A metal center is coordinated to **two atoms from the same amino acid residue**
* One of the coordinating atoms belongs to the **amino acid backbone**

This update improves the robustness of residue recognition, connectivity assignment, and bonded topology generation for complex coordination environments frequently observed in catalytic and structural metalloproteins.

---

## ⚡ Automated `tleap` Input Generation

The workflow is now fully automated for the generation of `tleap` input files.

### 🔧 Multiple Metal Centers

For metalloproteins containing **multiple metal centers** that are parameterized sequentially:

* Keep the generated `tleap.input` file for each individual metal complex
* After completing all parameterization steps, the corresponding `tleap` inputs can later be combined together into the final system preparation workflow

This provides greater flexibility for large multinuclear metalloproteins and mixed-metal systems.

---

## 🛠️ Bug Fixes & Stability Improvements

### Metalloprotein Workflow Fixes

Resolved several issues affecting metalloprotein parameterization, including:

* Improved handling of coordinated residues involving **backbone atoms**
* Correct processing of residues coordinating through **multiple atoms within the same amino acid**
* More reliable bonded topology generation for complex coordination geometries
* Enhanced consistency of automatically generated `tleap` inputs

### General Improvements

* Additional robustness improvements for automated workflow generation
* Improved handling of sequential metal parameterization procedures
* Minor stability and parsing fixes throughout the pipeline

---

## 🔄 Previously: easyPARM v4.20

**easyPARM v4.20** introduced major support for **metal–nucleic acid systems** and expanded metalloprotein force-field flexibility.

### 🧬 Metal–Nucleic Acid Support

easyPARM now fully supports:

* Metal-bound DNA and RNA
* Inner-sphere metal coordination to nucleobases and phosphate groups
* Metal-mediated base pairs and metallo-DNA architectures
* Nucleic acid–metal complexes and directed metal assemblies

### 🔧 Supported Nucleic Acid Force Fields

Users can select among:

* **OL15**
* **OL21**
* **OL24**
* **BSC0**
* **BSC1**
* **OL3**

📖 Tutorial:
[https://abdelazim-abdelgawwad.github.io/Tutorial/Tutorial/metallonucleicacid](https://abdelazim-abdelgawwad.github.io/Tutorial/Tutorial/metallonucleicacid)

### 🧪 Expanded Metalloprotein Force Fields

easyPARM now supports:

* **ff12SB**
* **ff14SB**
* **ff19SB**
* **AMBER-fb15**

### 🛠️ v4.20 Bug Fixes

* Fixed incorrect bond information in `Bond_info.dat`
* Improved handling of terminal residues coordinated to metals
* Fixed issues in MOL2 file generation
* Improved robustness of intermediate parameter files

---

## 🔄 Previously: easyPARM v4.15

### 🛠️ Bug Fixes & Debugging

* Fixed issues in Gaussian and ORCA output parsing
* Resolved problems with atom-type assignments
* Improved debug messages for clearer user feedback

### 📖 Documentation

* Added a dedicated metalloprotein parameterization tutorial

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

easyPARM integrates the RESP implementation by Alenaizan within the Psi4 framework.

### ⚡ Enhanced ORCA Support

* Native RESP fitting via `!RESP` in ORCA 6.0
* Optional easyPARM RESP workflow using user-defined ESP points
* Support for both ORCA-native and easyPARM-based RESP fitting

---

## 🧠 Previously: easyPARM v4.00

* Improved atom-type detection for challenging coordination environments
* Accurate recognition based on topology, hybridization, aromaticity, ring strain, and metal–ligand bonding
* Introduction of ff19SB typing for metal-coordinated standard residues
* Updated non-bonded metal parameters compatible with AMBER

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
> DOI: 10.1021/acs.jctc.4c01272

> Abdelazim M. A. Abdelgawwad and Antonio Francés-Monerris.
> *easyPARM v4.00: A Python-based tool for the automated parameterization of metalloproteins and metal–organic polyhedra with multiple metal centers.*
> *Journal of Chemical Physics*, 2025, **163**(22), 222501.
> DOI: 10.1063/5.0301038
