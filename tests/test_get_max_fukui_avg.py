import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Highlighting import get_max_fukui_avg

def test_get_max_fukui_avg():
    fukui_avg = {0: 0.1, 1: 0.3, 2: 0.2}
    max_idx, max_value = get_max_fukui_avg(fukui_avg)
    assert max_idx == 1, "Test failed: Incorrect max index"
    assert max_value == pytest.approx(0.3, rel=1e-9), "Test failed: Incorrect max value"

