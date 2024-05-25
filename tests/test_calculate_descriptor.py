import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from src.Rankchem.for_streamlit.Ranking import calculate_descriptor

def test_calculate_descriptor():
    smiles_list = ["CC=O", "CC(=O)C"]
    descriptor_type = 'N'
    iterations = 2
    descriptors = calculate_descriptor(smiles_list, descriptor_type, iterations)
    assert isinstance(descriptors, list), "Test failed: Result is not a list"
    assert len(descriptors) == len(smiles_list), "Test failed: Incorrect number of descriptors"

    for smiles, descriptor in descriptors:
        assert isinstance(smiles, str), "Test failed: SMILES is not a string"
        assert isinstance(descriptor, float), "Test failed: Descriptor is not a float"
