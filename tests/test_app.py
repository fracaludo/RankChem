import sys
import os
import pytest
from unittest.mock import MagicMock, patch

# Add the src directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

# Mock streamlit.session_state
@pytest.fixture(autouse=True)
def mock_streamlit_session_state(monkeypatch):
    # Create a mock for st.session_state
    mock_session_state = MagicMock()
    mock_session_state.page = "Home"
    
    # Patch the session_state attribute in the streamlit module
    monkeypatch.setattr("streamlit.session_state", mock_session_state)

# Import the functions after patching to ensure the patching works
from Rankchem.app import home, documentation, github, run_active_sites, run_ranking

def test_home():
    assert home is not None

def test_documentation():
    assert documentation is not None

def test_github():
    assert github is not None

def test_run_active_sites():
    assert run_active_sites is not None

def test_run_ranking():
    assert run_ranking is not None
