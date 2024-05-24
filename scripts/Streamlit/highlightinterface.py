import streamlit as st
from RankChem.Highlighting import smiles_to_xyz, calculate_fukui, average, get_max_fukui_avg


def visualize_molecule(smiles, max_idx, highlight_color, representation):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, rdDistGeom.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=800, height=600)
    view.addModel(mol_block, "mol")
    view.setStyle({representation.lower(): {}})
    view.addStyle({'model': -1, 'serial': max_idx-1}, {'sphere': {'radius': 0.5, 'color': highlight_color, 'opacity': 0.8}})
    view.zoomTo()
    showmol(view, height=500, width=500)

def main(smiles, fukui_type, num_iterations, highlight_color, representation):
    xyz_filename = smiles_to_xyz(smiles)
    elements, coordinates = read_xyz(xyz_filename)
    
    fukui_type_map = {"E": "local_nucleophilicity", "N": "local_electrophilicity"}
    fukui_type_str = fukui_type_map.get(fukui_type.upper())
    
    fukui_dicts = []
    for _ in range(num_iterations):
        fukui_dict = calculate_fukui(elements, coordinates, fukui_type_str)
        fukui_dicts.append(fukui_dict)
    
    fukui_avg = average(fukui_dicts)
    max_idx, max_value = get_max_fukui_avg(fukui_avg)
    
    st.write("The results are shown below:")
    st.write(f"The Fukui average values for each atom: {fukui_avg}")
    st.write(f"The Atom with max average Fukui value: Atom {max_idx} with a value of {max_value}")
    visualize_molecule(smiles, max_idx, highlight_color, representation)



st.sidebar.title("Navigation")
if st.sidebar.button("Documentation"):
    selected_page = "Documentation"
elif st.sidebar.button("Github"):
    selected_page = "Github"
elif st.sidebar.button("Active sites"):
    selected_page = "Active sites"
else:
    selected_page = "Active sites"

if selected_page == "Documentation":
    st.title("Documentation")
    st.write("Here you can provide documentation and usage instructions for your application.")
    # Add more content for documentation here

elif selected_page == "Github":
    st.title("Github")
    st.write("Provide a link to your Github repository.")
    st.markdown("[Github Repository](https://github.com/fracaludo/RankChem.git)")
    # Add more content for Github here

elif selected_page == "Active sites":
    st.title("Molecule Visualization and Fukui Analysis")
    st.write("Use the options in the sidebar to visualize the molecule and highlight Fukui indices.")

    # Input options
    user_input = st.sidebar.text_input("Enter a SMILES string:", "CC=O")
    fukui_type = st.sidebar.radio("Highlight nucleophile sites (N) or electrophile sites (E)?", ["E", "N"])
    num_iterations = st.sidebar.slider("Number of iterations", min_value=100, max_value=400, value=1)
    highlight_color = st.sidebar.color_picker("Highlight color", "#FF5733")
    representation = st.sidebar.selectbox("Choose Py3Dmol representation", ["Stick", "Line", "Cross", "Sphere", "Cartoon"])
    
    if st.sidebar.button("Submit"):
        main(user_input, fukui_type, num_iterations, highlight_color, representation)