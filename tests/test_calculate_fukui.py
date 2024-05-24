import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.Highlighting import calculate_fukui

def test_calculate_fukui():
    elements = ["C", "O", "H", "H"]
    coordinates = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.2], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
    fukui_type = "nucleophilicity"
    fukui_values = calculate_fukui(elements, coordinates, fukui_type)
    assert isinstance(fukui_values, dict), "Test failed: Result is not a dictionary"
    assert all(isinstance(key, int) and isinstance(value, float) for key, value in fukui_values.items()), "Test failed: Dictionary keys are not integers or values are not floats"
