import os
from rdkit import Chem
from rdkit.Chem import AllChem
from morfeus import read_xyz, XTB
import streamlit as st

def smiles_to_xyz(smiles, filename):
    '''
    Converts a given SMILES representation of a molecule into an XYZ file format.
    It uses RDKit to generate the 3D coordinates of the molecule and writes them into an XYZ file.

    Inputs:
    smiles: The SMILES representation of the molecule.
    filename: The name of the output XYZ file.

    Output:
    None
    '''
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    with open(filename, 'w') as f:
        f.write(Chem.MolToXYZBlock(mol))

def calculate_descriptor(smiles_list, descriptor_type, iterations):
    '''
    Calculates the specified global descriptor (nucleophilicity or electrophilicity) for a list of molecules.
    The calculation is performed multiple times for each molecule to obtain an average value.

    Inputs:
    smiles_list: A list of SMILES strings representing the molecules.
    descriptor_type: The type of descriptor to calculate ('N' for nucleophilicity, 'E' for electrophilicity).
    iterations: The number of iterations to perform for averaging the descriptor values.

    Output:
    Returns a list of tuples, where each tuple contains a SMILES string and its corresponding average descriptor value.
    '''
    descriptors = []

    for smiles in smiles_list:
        xyz_file = f"{smiles}.xyz"
        smiles_to_xyz(smiles, xyz_file)
        
        elements, coordinates = read_xyz(xyz_file)
        xtb = XTB(elements, coordinates)

        descriptor_values = []
        for _ in range(iterations):
            if descriptor_type == 'N':
                descriptor = xtb.get_global_descriptor("nucleophilicity", corrected=True)
            elif descriptor_type == 'E':
                descriptor = xtb.get_global_descriptor("electrophilicity", corrected=True)
            descriptor_values.append(descriptor)
        
        # Average the descriptor values over the iterations
        avg_descriptor = sum(descriptor_values) / len(descriptor_values)
        descriptors.append((smiles, avg_descriptor))
        
        # Clean up the temporary XYZ file
        os.remove(xyz_file)
    
    return descriptors

def rank_descriptors(descriptors):
    '''
    Ranks a list of molecules based on their descriptor values in descending order.

    Input:
    descriptors: A list of tuples, where each tuple contains a SMILES string and its corresponding descriptor value.

    Output:
    Returns a list of tuples, where each tuple contains a SMILES string and its corresponding descriptor value, sorted in descending order.
    '''
    ranked_descriptors = sorted(descriptors, key=lambda x: x[1], reverse=True)
    return ranked_descriptors

def run_Ranking():
    '''
    Main function for the Streamlit application.
    Provides a user interface for inputting a list of SMILES strings, selecting a descriptor type, and running the ranking process.

    Inputs:
    None

    Output:
    Displays the ranked list of molecules based on their descriptor values.
    '''
    st.title("Molecule Ranking")
    st.write("Use the options below to rank molecules in order of decreasing reacivity (Top molecule will be most reactive and bottom molecule the least reactive).")

    # Input widgets for Streamlit
    user_input = st.text_input("Enter a list of SMILES strings separated by commas:", "CC=O,CC(=O)C,O=COC,CN(C)C(C)=O")
    descriptor_type = st.radio("Should we be ranking according to nucleophilicity (N) or electrophilicity (E)?", ('N', 'E'))
    iterations = st.slider("Number of iterations (increase for a more precise result).", min_value=1, max_value=100, value=10)
    
    if st.button("Run Ranking"):
        smiles_list = [smiles.strip() for smiles in user_input.split(',')]
        descriptors = calculate_descriptor(smiles_list, descriptor_type, iterations)
        ranked_descriptors = rank_descriptors(descriptors)

        st.write("Ranked list of molecules based on their descriptors:")
        for smiles, descriptor in ranked_descriptors:
            st.write(f"{smiles}: {descriptor}")

if __name__ == "__main__":
    run_ranking()
