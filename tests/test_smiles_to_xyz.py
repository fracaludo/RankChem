import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Highlighting import smiles_to_xyz
import os

def test_smiles_to_xyz():
    filename = "test_molecule.xyz"
    result = smiles_to_xyz("C(=O)O", filename)
    assert result == filename, "Test failed: Incorrect filename returned"
    assert os.path.exists(filename), "Test failed: File not created"
    with open(filename, "r") as f:
        content = f.readlines()
    assert len(content) > 2, "Test failed: File content is too short"
    os.remove(filename)


