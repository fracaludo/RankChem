import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Ranking import calculate_descriptor

def test_calculate_descriptor():
    smiles_list = ["CC=O", "CC(=O)C"]
    descriptor_type = "N"
    iterations = 10
    descriptors = calculate_descriptor(smiles_list, descriptor_type, iterations)
    assert isinstance(descriptors, list), "Test failed: Result is not a list"
    assert len(descriptors) == len(smiles_list), "Test failed: Incorrect number of descriptors calculated"

