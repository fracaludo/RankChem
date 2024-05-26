import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from morfeus import read_xyz, XTB
from stmol import showmol
import py3Dmol

def smiles_to_xyz(smiles, filename="molecule.xyz"):
    '''
    This function converts a given SMILES representation of a molecule into its corresponding XYZ file format.
    It uses the RDKit library to generate the 3D coordinates of the molecule and writes them into an XYZ file.

    Inputs:
    smiles: The SMILES representation of the molecule.
    filename: (Optional) The name of the output XYZ file. Default is "molecule.xyz".

    Output:
    Returns the filename of the generated XYZ file.
    '''
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, rdDistGeom.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    with open(filename, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
    return filename

def calculate_fukui(elements, coordinates, fukui_type):
    '''
    This function computes the Fukui function values for each atom in a molecule.
    It utilizes quantum chemical calculations to determine the reactivity descriptors based on the molecule's electronic structure.

    Inputs:
    elements: A list of atomic symbols for the atoms in the molecule.
    coordinates: A list of lists containing the Cartesian coordinates (x, y, z) for each atom.
    fukui_type: The type of Fukui function to be calculated (e.g., nucleophilicity or electrophilicity).

    Output:
    Returns a dictionary mapping each atom index to its corresponding Fukui function value.
    '''
    xtb = XTB(elements, coordinates)
    fukui_values = xtb.get_fukui(fukui_type)
    return fukui_values

def average(fukui_dicts):
    '''
    This function computes the average Fukui function values across multiple dictionaries,
    each containing Fukui function values for the same set of atoms but for different molecular configurations or environments.

    Input:
    fukui_dicts: A list of dictionaries, where each dictionary maps atom indices to Fukui function values.

    Output:
    Returns a dictionary containing the average Fukui function values for each atom index.
    '''
    avg_fukui = {}
    num_dicts = len(fukui_dicts)
    
    for fukui_dict in fukui_dicts:
        for atom_idx, value in fukui_dict.items():
            if atom_idx not in avg_fukui:
                avg_fukui[atom_idx] = 0.0
            avg_fukui[atom_idx] += value
    
    for atom_idx in avg_fukui:
        avg_fukui[atom_idx] /= num_dicts
    
    return avg_fukui

def get_max_fukui_avg(fukui_avg):
    '''
    This function identifies the atom with the highest average Fukui function value from a given dictionary of average Fukui function values.

    Input:
    fukui_avg: A dictionary containing the average Fukui function values for each atom index.

    Output:
    Returns the index of the atom with the highest average Fukui function value and the corresponding value.
    '''
    max_idx = max(fukui_avg, key=fukui_avg.get)
    max_value = fukui_avg[max_idx]
    return max_idx, max_value

def visualize_molecule(smiles, max_idx, highlight_color, representation):
    '''
    This function visualizes a molecule using 3D molecular visualization tools.
    It utilizes RDKit to generate a molecular model and Py3Dmol to render the model in a 3D viewer.

    Inputs:
    smiles: The SMILES representation of the molecule.
    max_idx: The index of the atom with the highest average Fukui function value, used for visualization highlighting.
    highlight_color: The color to highlight the atom with the highest Fukui value.
    representation: The representation style for the 3D visualization (e.g., Stick, Line, Cross, Sphere).

    Output:
    Displays the interactive 3D visualization of the molecule.
    '''
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, rdDistGeom.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=1000, height=800) 
    view.addModel(mol_block, "mol")
    view.setStyle({representation.lower(): {}})
    view.addStyle({'model': -1, 'serial': max_idx-1}, {'sphere': {'radius': 0.5, 'color': highlight_color, 'opacity': 0.8}})
    view.zoomTo()
    showmol(view, height=600, width=800)  


def run_Highlighting():
    '''
    This function serves as the main entry point for the Streamlit application.
    It provides a user interface for entering a SMILES string, selecting options for Fukui analysis,
    and visualizing the molecule with highlighted active sites.

    Inputs:
    None

    Output:
    Displays the user interface and results of the Fukui analysis and molecular visualization.
    '''
    st.title("Molecule Visualization and Fukui Analysis")
    st.write("Use the options below to visualize the molecule and highlight Fukui indices.")

    # Input options moved from sidebar to main content
    user_input = st.text_input("Enter a SMILES string that you want to analyse:", "CC=O")
    fukui_type = st.radio("Highlight nucleophile sites (N) or electrophile sites (E)?", ["E", "N"])
    num_iterations = st.slider("Number of iterations (increase for a more precise result).", min_value=1, max_value=10, value=1)
    highlight_color = st.color_picker("Highlight color", "#FF5733")
    representation = st.selectbox("Choose a Py3Dmol representation style", ["Stick", "Line", "Cross", "Sphere"])
    
    if st.button("Submit"):
        # Generate XYZ file
        xyz_filename = smiles_to_xyz(user_input)
        elements, coordinates = read_xyz(xyz_filename)
        
        # Calculate Fukui indices
        fukui_type_map = {"E": "local_nucleophilicity", "N": "local_electrophilicity"}
        fukui_type_str = fukui_type_map.get(fukui_type.upper())
        
        fukui_dicts = []
        for _ in range(num_iterations):
            fukui_dict = calculate_fukui(elements, coordinates, fukui_type_str)
            fukui_dicts.append(fukui_dict)
        
        fukui_avg = average(fukui_dicts)
        max_idx, max_value = get_max_fukui_avg(fukui_avg)
        
        st.write(f"The Fukui average values are : {fukui_avg}")
        st.write(f"The atom with the highest average Fukui value is: {max_idx} with a value of {max_value}")
        
        visualize_molecule(user_input, max_idx, highlight_color, representation)
