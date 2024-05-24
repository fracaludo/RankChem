import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.Highlighting import main

def test_main():
    smiles = "C(=O)O"
    fukui_type = "E"
    num_iterations = 1
    main(smiles, fukui_type, num_iterations)
