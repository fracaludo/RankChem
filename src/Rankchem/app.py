
import streamlit as st
from pathlib import Path
import os
import sys

# Set up the path to import the scripts from the "for streamlit" directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'for streamlit')))

from Ranking import run_Ranking
from Highlighting import run_Highlighting 

# Paths to images
path = Path(__file__).parent / 'images'
path_img1 = path / 'rank.png'
path_img2 = path / 'rank.png'
path_img3 = path / 'rank.png'
path_img4 = path / 'rank.png'

# Initialize session state
if 'page' not in st.session_state:
    st.session_state.page = "Home"

def home():
    '''
    Renders the home page of the Streamlit application.
    Displays images and provides navigation links to other sections.

    Inputs:
    None

    Output:
    Displays the home page content.
    '''
    st.title("Home Page")
    st.write("Welcome to the Home Page. Use the images or the sidebar to navigate to different sections.")
    
    col1, col2 = st.columns(2)

    with col1:
        st.image(str(path_img1), use_column_width=True)
        if st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Documentation" target="_self">Documentation</a></h3>''', unsafe_allow_html=True):
            st.query_params["page"] = "Documentation"

        st.markdown("<p style='text-align: center;'>This section provides documentation and usage instructions.<br />&nbsp;</p>", unsafe_allow_html=True)

        st.image(str(path_img3), use_column_width=True)
        if st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Highlighting" target="_self">Highlighting</a></h3>''', unsafe_allow_html=True):
            st.query_params["page"] = "Highlighting"

        st.markdown("<p style='text-align: center;'>Visualize molecule and highlight active sites.<br />&nbsp;</p>", unsafe_allow_html=True)

    with col2:
        st.image(str(path_img4), use_column_width=True)
        if st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Ranking" target="_self">Ranking</a></h3>''', unsafe_allow_html=True):
            st.query_params["page"] = "Ranking"

        st.markdown("<p style='text-align: center;'>Rank molecules based on descriptors.<br />&nbsp;</p>", unsafe_allow_html=True)

def documentation():
    '''
    Renders the documentation page of the Streamlit application.
    Provides documentation and usage instructions for the application.

    Inputs:
    None

    Output:
    Displays the documentation content.
    '''
    st.title("Documentation")
    st.write("Here you can provide documentation and usage instructions for your application.")



# Sidebar navigation with buttons
st.sidebar.title("RankChem")
if st.sidebar.button("Home"):
    st.session_state.page = "Home"
    st.query_params["page"] = "Home"
if st.sidebar.button("Documentation"):
    st.session_state.page = "Documentation"
    st.query_params["page"] = "Documentation"
if st.sidebar.button("Highlighting"):
    st.session_state.page = "Highlighting"
    st.query_params["page"] = "Highlighting"
if st.sidebar.button("Ranking"):
    st.session_state.page = "Ranking"
    st.query_params["page"] = "Ranking"


# Add GitHub badge link at the bottom of the sidebar
st.sidebar.markdown(
    '''
    ---
    [![jhc github](https://img.shields.io/badge/GitHub-181717.svg?style=flat&logo=github)](https://github.com/fracaludo/RankChem.git)
    '''
)

# Navigation based on session state
if st.session_state.page == "Home":
    home()
elif st.session_state.page == "Documentation":
    documentation()
elif st.session_state.page == "Highlighting":
    run_Highlighting()
elif st.session_state.page == "Ranking":
    run_Ranking()

