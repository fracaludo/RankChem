import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from src.Rankchem.app import home, documentation, github, run_active_sites, run_ranking

def test_home():
    home()
    # Assuming no return value to assert

def test_documentation():
    documentation()
    # Assuming no return value to assert

def test_github():
    github()
    # Assuming no return value to assert

def test_run_active_sites():
    run_active_sites()
    # Assuming no return value to assert

def test_run_ranking():
    run_ranking()
    # Assuming no return value to assert

if __name__ == "__main__":
    pytest.main()
