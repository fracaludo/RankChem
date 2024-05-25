import sys
import os
import pytest

# Add the src directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

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
