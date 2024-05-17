# ChemInterface package for reactivity analysis
### Project in practical programming in chemistry course EPFL CH-200

## 🔥 Usage

This repository provides the user with a package which will display an interactive interface where the user will be able to input the SMILEs of multiple molecules. The user will then be able choose if he wants to analyze their electrophilicity or nucleophilicity. The interfcace will then display the reactivity ranking of the molecules based on the option you chose. Moreover, it will display their 3D structure with their most nucleophilic and electrophilic sites highlighted.

Developpers:
- Ludovica Fracassi, https://github.com/fracaludo
- Emma Kappeler, https://github.com/kappelemma


Before installing everything, let's firts define electrophilicity and nucleophilicity!!

**What is electrophilicity and nucleophilicity ?**
   - Electrophilicity and nucleophilicity are mesures of the reactivity of molecules....

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


