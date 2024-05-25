import streamlit as st
from pathlib import Path
import os
import sys

# Set up the path to import the scripts from the "for streamlit" directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'for_streamlit')))

# Import the run_ranking function from ranking.py
from Ranking import run_Ranking

# Import the run_active_sites function from active_sites.py
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
        st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Documentation" target="_self">Documentation</a></h3>''', unsafe_allow_html=True)
        st.markdown("<p style='text-align: center;'>This section provides documentation and usage instructions.<br />&nbsp;</p>", unsafe_allow_html=True)

        st.image(str(path_img3), use_column_width=True)
        st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Active Sites" target="_self">Active Sites</a></h3>''', unsafe_allow_html=True)
        st.markdown("<p style='text-align: center;'>Visualize molecule and highlight active sites.<br />&nbsp;</p>", unsafe_allow_html=True)

    with col2:
        st.image(str(path_img2), use_column_width=True)
        st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Github" target="_self">Github</a></h3>''', unsafe_allow_html=True)
        st.markdown("<p style='text-align: center;'>Link to the Github repository.<br />&nbsp;</p>", unsafe_allow_html=True)

        st.image(str(path_img4), use_column_width=True)
        st.markdown(f'''<h3 style='text-align: center;'><a href="?page=Ranking" target="_self">Ranking</a></h3>''', unsafe_allow_html=True)
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

def github():
    '''
    Renders the GitHub page of the Streamlit application.
    Provides a link to the GitHub repository of the application.

    Inputs:
    None

    Output:
    Displays the GitHub link.
    '''
    st.title("Github")
    st.write("This is RankChem's GitHub repo:)")
    st.markdown("[Github Repository](https://github.com/fracaludo/RankChem.git)")

# Check URL parameters for page navigation
query_params = st.experimental_get_query_params()
if 'page' in query_params:
    page_from_query = query_params['page'][0]
    if page_from_query in ["Home", "Documentation", "Github", "Active Sites", "Ranking"]:
        st.session_state.page = page_from_query
    else:
        st.session_state.page = "Home"  # default to Home if the page is not valid

# Sidebar navigation with buttons
st.sidebar.title("RankChem")
if st.sidebar.button("Home"):
    st.session_state.page = "Home"
    st.experimental_set_query_params(page="Home")
if st.sidebar.button("Documentation"):
    st.session_state.page = "Documentation"
    st.experimental_set_query_params(page="Documentation")
if st.sidebar.button("Github"):
    st.session_state.page = "Github"
    st.experimental_set_query_params(page="Github")
if st.sidebar.button("Active Sites"):
    st.session_state.page = "Active Sites"
    st.experimental_set_query_params(page="Active Sites")
if st.sidebar.button("Ranking"):
    st.session_state.page = "Ranking"
    st.experimental_set_query_params(page="Ranking")

# Navigation based on session state
if st.session_state.page == "Home":
    home()
elif st.session_state.page == "Documentation":
    documentation()
elif st.session_state.page == "Github":
    github()
elif st.session_state.page == "Active Sites":
    run_Highlighting()
elif st.session_state.page == "Ranking":
    run_Ranking()
