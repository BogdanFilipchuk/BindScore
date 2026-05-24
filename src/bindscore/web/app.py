import streamlit as st
from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_data as fetch_pdb
import py3Dmol
import streamlit.components.v1 as components # if using regular st.html, gives error
import numpy as np
import matplotlib.colors as mcolors
#from import get_chains # to get the interaction dictionary
from bindscore.scoring.total_entropy import * # to get the entropy breakdown

# auxiliary functions and tools
interaction_color = {                               # Define colors for interaction types
        "hydrogen_bond": "#bb2727",
        "pi_pi_stacking": "#06e7c9",
        "disulfide_bond": "#ceca00",
        "metal_coordination": "#ec7d8c",
        "salt_bridge": "#762a83",
        "hydrophobic_contact": "#1d29d4",
        "halogen_bond": "#07d600",
        "dipole-dipole": "#036400",
        "unclassified": "#999999"}


class Atom(dict):
    def __init__(self, atom_dict):
        self["type"] = atom_dict["type"]
        self["atom_seq"] = atom_dict["atom_seq"]


    def __str__(self):
        line = list(" " * 80)

        line[0:6] = self["type"].ljust(6)
        line[6:11] = self["idx"].ljust(5)
        line[12:16] = self["name"].ljust(4)
        line[17:20] = self["resname"].ljust(3)
        line[22:26] = str(self["resid"]).ljust(4)
        line[30:38] = str(self["x"]).rjust(8)
        line[38:46] = str(self["y"]).rjust(8)
        line[46:54] = str(self["z"]).rjust(8)
        line[76:78] = self["sym"].rjust(2)
        return "".join(line) + "\n"
    
def get_interaction_color(interaction_type, mode = ["by_interaction", "by_binding_energy"], dH = None, norm = None):
    if mode == "by_binding_energy":                         # If coloring by binding energy differences
        cmap = mcolors.LinearSegmentedColormap.from_list(   # Create a custom colormap for the binding energy differences
            [(0.0, "#2166ac"), (0.5, "#a1a1a1"), (1.0, "#d6604d")],
            N=512,
        )
        
    
    return interaction_color.get(interaction_type, "#999999")  # Default to gray if type is unknown  
#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####TEMPERATURE SLIDER#####
#st.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### protein input and visualization
protein_input = fetch_pdb(st.text_input("Enter protein PDB RCSB database ID (e.g. 1A2B)", "1A2B"))  # input the protein PDB ID
# residues_importance = calculate_residue_importance(protein_input)

### selectbox with the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["cartoon", "surface", "stick"])  # the selectbox

### visualization button and viewer
if st.button("Visualize"):
    st.spinner(text="Loading...", width="content")
    viewer = py3Dmol.view(width=600, height=600)        # Create viewer
    viewer.addModel(protein_input, "pdb")               # PDB format
    if Protein_drawmode == "surface":                   # Surface drawing requires a separate call
        viewer.addSurface(py3Dmol.VDW, {"opacity": 1})
    elif Protein_drawmode == "cartoon":                   # Apply selected style
        viewer.setStyle({'cartoon': {'color': 'spectrum'}}) # Style options
    else:                                               # Default to stick style
        viewer.setStyle({'stick': {}})

### interaction calculation and gradient coloring  
if st.button("Calculate"):                              # button to load the protein viewer
    st.spinner(text="Calculating interactions...", width="content")      # Show spinner while processing
    interactions_full = func(protein_input)             # Call the function to calculate interactions and binding energy differences
    norm = mcolors.TwoSlopeNorm(                        # Normalize the binding energy around zero
            vmin=min(dH), vcenter=0.0, vmax=max(dH)
        )
        
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

### Interactions table 
# data = imported_interactions
# st.dataframe(data)  # Display interactions in a table format

### Entropy table
results = compute_total_entropy(
    complex_pdb = protein_input, 
    chain_a = chain_1, 
    chain_b = chain_2, 
    return_breakdown = True
)
results_dict = {"Entropy type": "Free Energy contribution (J/(mol))",
    "Translational-Rotational": results.dS_trans_rot, 
    "Hydrophobic": results.dS_hydrophobic, 
    "Side-chain": results.dS_sidechain, 
    "Backbone": results.dS_backbone, 
    "Total": results.dS_total
}
st.table(results_dict)