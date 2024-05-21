# ChemInterface package for reactivity analysis
### Project in practical programming in chemistry course EPFL CH-200
[![alt text](https://img.shields.io/badge/Python-14354C?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org)
![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)
[![Emma1](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)
![license](https://custom-icon-badges.demolab.com/github/license/denvercoder1/custom-icon-badges?logo=law&logoColor=white)
![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)
<img src="https://user-images.githubusercontent.com/7164864/217935870-c0bc60a3-6fc0-4047-b011-7b4c59488c91.png" alt="Streamlit logo" style="margin-top:20px"></img>

# 🔥 Usage

This repository provides the user with a package which will display an interactive interface where the user will be able to input the SMILEs of multiple molecules. The user will then be able choose if he wants to analyze their electrophilicity or nucleophilicity. The interfcace will then display the reactivity ranking of the molecules based on the option you chose. Moreover, it will display their 3D structure with their most nucleophilic and electrophilic sites highlighted.

Developpers:
- Ludovica Fracassi, https://github.com/fracaludo
- Emma Kappeler, https://github.com/kappelemma


Before installing everything, let's first define electrophilicity and nucleophilicity!!

**What is electrophilicity and nucleophilicity ?**
   - Electrophilicity and nucleophilicity are measures of the reactivity of molecules....

Now let's go through the steps required to use this package.

## 👩‍💻 Installation

RankChem can be installed using pip as followed:
```bash
pip install rankchem
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
```
Specifically, from these packages, the following subpackages are required:
```bash
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from morfeus import read_xyz, XTB
import py3Dmol
```

## 🎥 How it works

....

## 🚥 Let's get started!

....
