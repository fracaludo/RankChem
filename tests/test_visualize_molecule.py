import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from Rankchem.for_streamlit.Highlighting import visualize_molecule

def test_visualize_molecule():
    smiles = "C(=O)O"
    adj_max_idx = 1
    highlight_color = "red"  
    representation = "stick"  
    visualize_molecule(smiles, adj_max_idx, highlight_color, representation)



