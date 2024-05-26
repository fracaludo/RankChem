![logo](nkChem.png)

# ChemInterface package for reactivity analysis
### Project in practical programming in chemistry course EPFL CH-200
[![alt text](https://img.shields.io/badge/Python-14354C?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org)
![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)
[![Emma1](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)
![license](https://img.shields.io/badge/License-MIT-ac8b11.svg?style=for-the-badge&labelColor=blue)

# 🔥 Usage

This repository provides the user with a package which will display an interactive interface where the user will be able to input the SMILEs of multiple molecules. The user will then be able choose if he wants to analyze their electrophilicity or nucleophilicity. The interfcace will then display the reactivity ranking of the molecules based on the option you chose. Moreover, it will display their 3D structure with their most nucleophilic and electrophilic sites highlighted.

The developpers of this package are:
- Ludovica Fracassi [![jhc github](https://img.shields.io/badge/GitHub-fracaludo-181717.svg?style=flat&logo=github)](https://github.com/fracaludo)
- Emma Kappeler [![jhc github](https://img.shields.io/badge/GitHub-kappelemma-181717.svg?style=flat&logo=github)](https://github.com/kappelemma)


Before installing everything, let's first define electrophilicity and nucleophilicity!!

**What is electrophilicity and nucleophilicity ?**

Electrophilicity and nucleophilicity are fundamental concepts in organic chemistry that describe the reactivity of molecules. Electrophilicity refers to the ability of a molecule or ion to accept an electron pair, making it an electron-loving species (electrophile). Electrophiles typically have a positive charge, partial positive charge, or an electron-deficient atom, making them attracted to electron-rich regions. On the other hand, nucleophilicity describes the ability of a molecule or ion to donate an electron pair, making it an electron-rich species (nucleophile). Nucleophiles are usually negatively charged or have lone pairs of electrons, such as anions, amines, and alcohols. The interaction between nucleophiles and electrophiles drives many chemical reactions, particularly in organic synthesis, where nucleophiles attack electrophiles to form new bonds.

Now let's go through the steps required to use this package.

## 👩‍💻 Installation

RankChem can be installed using pip as followed:
```bash
pip install Rankchem
```

The package can also be installed from source by running the following commands:

First, clone the repository from github and go in the folder.
```bash
git clone [https://github.com/fracaludo/PPchem-project.git]
cd rankchem
```
Then, install the package using :
```bash
pip install -e .
```
Or by following this single command:

```bash
pip install git+https://github.com/fracaludo/PPchem-project.git
```
#### Imported packages

In order to run the package correctly, the following packages need to be installed using the following commands. Moreover, this package works for python ....

```bash
conda config --add channels conda-forge

conda install -c conda-forge morfeus-ml
conda install "libblas=*=*mkl"
conda install xtb-python
conda install -c conda-forge rdkit
conda install -c conda-forge pyvistaqt
conda install -c conda-forge numpy
conda install -c conda-forge py3Dmol
conda install -c streamlit
conda install -c stmol


```
Specifically, from these packages, the following subpackages are required:
```bash
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from morfeus import read_xyz, XTB
import py3Dmol
import streamlit as st

```

## 🎥 How it works

The interface enables the user to input a smiles of a molecule and choose if nucleophilicity or electrophilicity is highlighted. The interface will then display the molecule in 3D with the chosen site highlighted. An example of the input and output is shown below:

XXX

## 🚥 Let's get started!

The streamlit interface can now finally be used!!

Just click on the logo below ;)

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://<your-custom-subdomain>.streamlit.app)
....
