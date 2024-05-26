import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Highlighting import calculate_fukui

def test_calculate_fukui():
    elements = ['C', 'O', 'H', 'H']
    coordinates = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    fukui_type = 'nucleophilicity'
    fukui_values = calculate_fukui(elements, coordinates, fukui_type)
    assert isinstance(fukui_values, dict), "Test failed: Result is not a dictionary"
    for key, value in fukui_values.items():
        assert isinstance(key, int), "Test failed: Key is not an integer"
        assert isinstance(value, float), "Test failed: Value is not a float"


