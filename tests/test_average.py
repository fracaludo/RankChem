import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from src.Rankchem.for_streamlit.Highlighting import average

def test_average():
    fukui_dicts = [
        {0: 0.1, 1: 0.2, 2: 0.3},
        {0: 0.2, 1: 0.3, 2: 0.4},
        {0: 0.3, 1: 0.4, 2: 0.5}
    ]
    avg_fukui = average(fukui_dicts)
    assert isinstance(avg_fukui, dict), "Test failed: Result is not a dictionary"
    assert avg_fukui[0] == pytest.approx(0.2, rel=1e-9), "Test failed: Incorrect average value for atom 0"
    assert avg_fukui[1] == pytest.approx(0.3, rel=1e-9), "Test failed: Incorrect average value for atom 1"
    assert avg_fukui[2] == pytest.approx(0.4, rel=1e-9), "Test failed: Incorrect average value for atom 2"

