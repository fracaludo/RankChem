import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Ranking import rank_descriptors

def test_rank_descriptors():
    descriptors = [('SMILES1', 0.5), ('SMILES2', 0.3), ('SMILES3', 0.7)]
    ranked_descriptors = rank_descriptors(descriptors)
    assert isinstance(ranked_descriptors, list), "Test failed: Result is not a list"
    assert len(ranked_descriptors) == len(descriptors), "Test failed: Number of elements in result is incorrect"
    assert ranked_descriptors == [('SMILES3', 0.7), ('SMILES1', 0.5), ('SMILES2', 0.3)], "Test failed: Incorrect ranking"
