import streamlit as st
from pathlib import Path
import os
import sys

st.set_page_config(layout="wide")

current_dir = Path(__file__).parent.resolve()
for_streamlit_dir = current_dir / 'for_streamlit'
sys.path.insert(0, str(for_streamlit_dir))


if not (for_streamlit_dir / 'Ranking.py').exists():
    st.error("Ranking.py does not exist in the 'for streamlit' directory.")
if not (for_streamlit_dir / 'Highlighting.py').exists():
    st.error("Highlighting.py does not exist in the 'for streamlit' directory.")

try:
    from Ranking import run_Ranking
    from Highlighting import run_Highlighting 
except ModuleNotFoundError as e:
    st.error(f"Error importing modules: {e}")
    run_Ranking = run_Highlighting = None



path_img2 = current_dir / 'images' / 'rank.png'
path_img1 = current_dir / 'images' / 'highlight.png'


if 'page' not in st.session_state:
    st.session_state.page = "Home"

def home():
    st.title("RankChem Homepage")
    st.write("Welcome to the RankChem Home Page. Use the buttons or the sidebar to navigate to the different sections :)")
    
    col1, col2 = st.columns(2)

    with col1:
        st.image(str(path_img1), use_column_width=True)
        if st.button("Go to Highlighting"):
            st.session_state.page = "Highlighting"

        st.markdown("<p style='text-align: center;'>Visualization of molecules and highlights of active sites.<br />&nbsp;</p>", unsafe_allow_html=True)

    with col2:
        st.image(str(path_img2), use_column_width=True)
        if st.button("Go to Ranking"):
            st.session_state.page = "Ranking"

        st.markdown("<p style='text-align: center;'>Ranking of molecules based on descriptors.<br />&nbsp;</p>", unsafe_allow_html=True)


# Sidebar navigation with buttons
st.sidebar.title("RankChem")
if st.sidebar.button("Home"):
    st.session_state.page = "Home"
if st.sidebar.button("Highlighting"):
    st.session_state.page = "Highlighting"
if st.sidebar.button("Ranking"):
    st.session_state.page = "Ranking"

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
elif st.session_state.page == "Highlighting":
    if run_Highlighting:
        run_Highlighting()
    else:
        st.error("Highlighting module not found.")
elif st.session_state.page == "Ranking":
    if run_Ranking:
        run_Ranking()
    else:
        st.error("Ranking module not found.")
