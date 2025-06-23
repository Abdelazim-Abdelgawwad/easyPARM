# 🚀 easyPARM v4.00 Released! 🎉

We’re excited to announce **easyPARM v4.00**, a major update packed with **critical enhancements** and **new force field features** for even more robust and intelligent parameter generation!

## 🔥 NEW in v4.00:

🧠 **Improved Atom Type Detection**
Our upgraded atom type assignment engine now features enhanced recognition of chemical environments, especially in **challenging coordination complexes**.
This ensures more **accurate atom typing** using **molecular topology**, **ring strain**, **hybridization**, and **aromaticity analysis**, with better **metal-ligand recognition**.

🧬 **ff19SB for Metal-Coordinated Standard Residues**
Standard amino acid residues **coordinated to metal centers** are now typed using **ff19SB**, improving structural realism and parameter consistency in **metal-binding proteins**.

🔬 **Updated Non-Bonded Parameters for metals**  
  includes updated **Lennard-Jones (van der Waals) parameters** for metal:
  - **Monovalent metal ions** from [Li et al., *J. Chem. Inf. Model.* (2021)](https://doi.org/10.1021/acs.jcim.0c01390)
  - **Divalent metal ions** from [Zhang et al., *J. Chem. Theory Comput.* (2021)](https://doi.org/10.1021/acs.jctc.0c00194)  
  These updates provide better consistency with **AMBER-compatible force fields** and improve the reliability of **ion–ligand interaction modeling**.

🚨 **Bug Fixes**
Numerous bugs have been fixed for better stability. 

---

## 🔄 Previously in v3.30:

🔍 **Enhanced GAFF Atom Type Assignment**
Advanced fallback logic for **GAFF atom typing** when **Antechamber fails**, particularly for **multi-metal and organometallic compounds**.

## ✅ Also in v3.25 & v3.20:

* **Non-Interactive Mode** for batch processing
* **Unlimited Multi-metal Support**
* **Metalloprotein Library Integration**

---

## 📚 Documentation

📖 [Tutorial & Documentation](https://abdelazim-abdelgawwad.github.io/Tutorial/)
📂 Explore usage examples, configuration tips, and more on the website!

---

## 💡 Contribute & Suggest Improvements

Have ideas to make **easyPARM** even better?
We welcome all **issues**, **pull requests**, and **community discussions**:

* [Open an Issue](https://github.com/Abdelazim-Abdelgawwad/easyPARM/issues)
* [Submit a Pull Request](https://github.com/Abdelazim-Abdelgawwad/easyPARM/pulls)
* Join the conversation in the **Discussions** tab

Let’s shape the future of molecular parameterization together! 🧪⚡

---

## 📖 How to Cite

If you use **easyPARM** in your research, please cite:

> Abdelazim M. A. Abdelgawwad and Antonio Francés-Monerris.
> *easyPARM: Automated, Versatile, and Reliable Force Field Parameters for Metal-Containing Molecules with Unique Labeling of Coordinating Atoms.*
> *J. Chem. Theory Comput.*, Article ASAP.
> [DOI: 10.1021/acs.jctc.4c01272](https://doi.org/10.1021/acs.jctc.4c01272)

---

