import streamlit as st
# import matplotlib.pyplot as plt
from bindscore.pdb_file_treatment.pdb_utils import fetch_pdb_data as fetch_pdb
import py3Dmol
import streamlit.components.v1 as components # if using regular st.html, gives error

#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####SIDEBAR#####
#st.sidebar.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### protein input and calculations
protein_input = fetch_pdb(st.text_input("Enter protein PDB RCSB database ID (e.g. 1A2B)", "1A2B"))  # input the protein PDB ID
# residues_importance = calculate_residue_importance(protein_input)

### selectbox with the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["cartoon", "surface", "stick"])  # the selectbox

### protein viewer
if st.button("Visualize"):                              # button to load the protein viewer
    st.spinner(text="Loading...", width="content")      # Show loading spinner while processing
    viewer = py3Dmol.view(width=600, height=600)        # Create viewer
    viewer.addModel(protein_input, "pdb")               # PDB format
    if Protein_drawmode == "surface":                   # Surface drawing requires a separate call
        viewer.addSurface(py3Dmol.VDW, {"opacity": 1})
    elif Protein_drawmode == "cartoon":                   # Apply selected style
        viewer.setStyle({'cartoon': {'color': 'spectrum'}}) # Style options
    else:                                               # Default to stick style
        viewer.setStyle({'stick': {}})                   # Stick style
    # for residue in residues_importance:                 # Highlight important residues (if calculated)
    #     color = residues_importance[residue]            # Get importance score for the residue
    #     if Protein_drawmode == "surface":               # Surface drawing requires a separate call
    #         viewer.addSurface(py3Dmol.VDW, {"opacity": 1}, color=color)
    #     elif Protein_drawmode == "cartoon":                   # Apply selected style
    #         viewer.setStyle({'resi': residue},  {'cartoon': {'color': color}}) # Style options
    viewer.zoomTo()                                     # Zoom to fit the molecule
    html = viewer._make_html()                          # Render in Streamlit
    components.html(html, width=600, height=600)        # Size of the streamlit component

### Interactions 
# data = imported_interactions
# st.dataframe(data)  # Display interactions in a table format

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