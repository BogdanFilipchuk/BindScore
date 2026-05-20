import streamlit as st
import matplotlib.pyplot as plt
<<<<<<< HEAD
=======
import numpy as np
>>>>>>> ed0c948053e0d6a53a645f320b92f2c597b8e464
from bindscore.pdb_file_treatment.pdb_utils import fetch_pdb_data as fetch_pdb
import py3Dmol
import streamlit.components.v1 as components

#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####SIDEBAR#####
#st.sidebar.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### protein input
protein_input = fetch_pdb(st.text_input("Enter protein PDB RCSB database ID (e.g. 1A2B)", "1A2B"))  # input the protein PDB ID

### selectbox with the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["stick", "cartoon", "surface"])  # the selectbox
py3dmol_drawmode = {"stick": {"stick": {}}, "cartoon": {"cartoon": {"color": "spectrum"}}} # dictionary to map the selectbox options to py3Dmol styles

### protein viewer
if st.button("Load Protein"):                           #button to load the protein viewer
    viewer = py3Dmol.view(width=600, height=600)        # Create viewer
    viewer.addModel(protein_input, "pdb")
    if Protein_drawmode == "surface":                   # Surface drawing requires a separate call
        viewer.addSurface(py3Dmol.VDW, {"opacity": 1})
    else:
        viewer.setStyle(py3dmol_drawmode[Protein_drawmode]) # Style options
    viewer.zoomTo()                                     # Zoom to fit the molecule
    html = viewer._make_html()                          # Render in Streamlit
    components.html(html, width=600, height=600)        # Size of the streamlit component

### pH Slider (expander)
# with st.expander("pH Scale"):
#     fig, ax = plt.subplots(figsize=(10, 0.2))           # Create a wide, short figure (horizontal line) for the pH scale
#     gradient = np.linspace(0, 1, 256).reshape(1, -1)    # Create gradient
#     ax.imshow(                                          #Show the gradient in graph
#         gradient,
#         aspect='auto',
#         cmap='RdYlBu',
#         extent=[0.0, 14.0, 0.0, 1.0]
#     )
#     ph = st.slider("Choose pH", 0.0, 14.0, 7.0)
#     ax.set_yticks([])                                   # Labels
#     ax.set_xticks(range(15)) 
#     ax.axvline(ph, color='black', linewidth=2)          # Indicator line
#     st.pyplot(fig)
#     if ph < 6:                                          # pH Interpretation
#         st.error("Acidic")
#     elif 6 <= ph < 7:
#         st.success("Neutral")
#     else:
#         st.info("Basic")