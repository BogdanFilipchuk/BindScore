import streamlit as st
from bindscore.pdb_file_treatment.pdb_utils_fetch import fetch_pdb_data as fetch_pdb
import py3Dmol
import streamlit.components.v1 as components # if using regular st.html, gives error
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
# from bindscore.scoring.total_entropy import * # to get the entropy breakdown
from bindscore.parsing.pdb_utils_protein import Protein_Structure # to get the chains
from bindscore.parsing.enthalpy_run import get_dataset_interaction_list # to get the interactions and binding energy differences

# auxiliary functions and tools
chain_colors = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231',
    '#911eb4', '#42d4f4', '#f032e6', '#bfef45',
    '#fabed4', '#469990', '#dcbeff', '#9A6324',
]
interaction_color = {                               # Define colors for interaction types
        "hydrogen_bond": "#bb2727",
        "pi_pi_stacking": "#06e7c9",
        "disulfide_bond": "#ceca00",
        "metal_coordination": "#ec7d8c",
        "salt_bridge": "#762a83",
        "hydrophobic_contact": "#1d29d4",
        "halogen_bond": "#07d600",
        "dipole-dipole": "#036400",
        "unclassified": "#999999",
}

def get_interaction_color(energy, norm):                       
    cmap = mcolors.LinearSegmentedColormap.from_list("bgr",  # Create a custom colormap for the binding energy differences
        [(0.0, "#2166ac"), (0.5, "#a1a1a1"), (1.0, "#d6604d")],
        N=512,
    )
    return mcolors.to_hex(cmap(norm(energy)))  # Map the energy to a color using the colormap and normalization

def get_binding_energy_summary(interactions: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(interactions)
    summary = (
        df.groupby("interaction type")["binding energy"]
        .sum()
        .reset_index()
        .rename(columns={"interaction type": "Interaction Type", "binding energy": "Total Binding Energy"})
        .sort_values("Total Binding Energy")
    )
    return summary

#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####TEMPERATURE SLIDER#####
#st.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### protein input and visualization
protein_input_name = st.text_input("Enter protein PDB RCSB database ID (e.g. 1A2B)", "1A2B")
protein_input = fetch_pdb(protein_input_name)  # input the protein PDB ID
protein_input_parsed = Protein_Structure(protein_input)  # turn data into a Protein_Structure object
chains = protein_input_parsed.chains()  # get the chains for interaction analysis
# residues_importance = calculate_residue_importance(protein_input)

### selectbox with the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["cartoon", "stick"])  # the selectbox

### visualization button and viewer
if st.button("Visualize"):
    st.spinner(text="Loading...", width="content")
    
    viewer = py3Dmol.view(width=750, height=750)        # Create viewer
    viewer.addModel(protein_input, "pdb")               # PDB format
    
    chains_with_colors = {chain: color for chain, color in zip(chains, chain_colors)}  # Map chains to colors

    # if Protein_drawmode == "surface":                 # Surface drawing requires a separate call, removed for simplicity 
    #     viewer.addSurface(py3Dmol.VDW, {"opacity": 1})
    
    for chain in chains:                                # Apply colors to each chain according to the visualization mode
        if Protein_drawmode == "cartoon":                                
            viewer.setStyle({'chain': chain}, {'cartoon': {'color': chains_with_colors[chain]}})
        else:                                               
            viewer.setStyle({'chain': chain}, {'stick': {'color': chains_with_colors[chain]}})
    
    viewer.zoomTo()                                     # Zoom to fit the molecule
    
    html = viewer._make_html()                          # Render in Streamlit
    components.html(html, width=750, height=750)        # Size of the streamlit component
    
    patches = [mpatches.Patch(color=color, label=f'Chain {chain}') for chain, color in chains_with_colors.items()] # Create legend patches for each chain
    fig, ax = plt.subplots(figsize=(1, 0.2 * len(chains_with_colors)))
    ax.legend(handles=patches, loc='center', fontsize=12)
    ax.axis('off')
    st.pyplot(fig)

### choosing the chains for interaction analysis
st.write("Select two chains for interaction analysis")  # Instruction for the user
col1,col2 =st.columns(2)                                # Create two columns for chain selection
with col1:
    chain_1 = st.selectbox("Select chain A for interaction analysis", chains)
with col2:
    chain_2 = st.selectbox("Select chain B for interaction analysis", chains)

### interaction calculation and gradient coloring  
if st.button("Analyze the interactions"):               # Button to start the interaction analysis
    st.spinner(text="Calculating...", width="content")  # Show spinner while processing
    interactions_full = get_dataset_interaction_list(protein_input_name, chain_1, chain_2)# Calculate interactions and binding energy differences
    energy_values = []
    for interaction in interactions_full:
        energy_values.append(interaction["binding energy"])
    norm = mcolors.TwoSlopeNorm(                        # Normalize the binding energy around zero
            vmin=min(energy_values), vcenter=0.0, vmax=max(energy_values))
    
    viewer = py3Dmol.view(width=750, height=750)        # Create viewer
    viewer.addModel(protein_input, "pdb")               # PDB format

    for interaction in interactions_full:                          
        if Protein_drawmode == "cartoon":                                
            viewer.setStyle( {'cartoon': {'color': "black"}}) # Reset all to white first to make the interaction colors stand out
            viewer.setStyle({'resi': interaction["atom1"]["residue_seq"]}, {'cartoon': {'color': get_interaction_color(interaction["binding energy"], norm) }})
            viewer.setStyle({'resi': interaction["atom2"]["residue_seq"]}, {'cartoon': {'color': get_interaction_color(interaction["binding energy"], norm) }})
        else:                                               
            viewer.setStyle( {'stick': {'color': "black"}}) # Reset all to white first to make the interaction colors stand out
            viewer.setStyle({'resi': interaction["atom1"]["residue_seq"]}, {'stick': {'color': get_interaction_color(interaction["binding energy"], norm) }})
            viewer.setStyle({'resi': interaction["atom2"]["residue_seq"]}, {'stick': {'color': get_interaction_color(interaction["binding energy"], norm) }})

    viewer.zoomTo()                                     # Zoom to fit the molecule
    html = viewer._make_html()                          # Render in Streamlit
    components.html(html, width=750, height=750)        # Size of the streamlit component

### Interactions table
    enthalpy_data = get_binding_energy_summary(interactions_full)  # Prepare the interaction data for display in a table format
    st.dataframe(enthalpy_data)                                    # Display interactions in a table format

### Entropy table - the entropy calculations turned out to be very complicated.
# results = compute_total_entropy(
#     complex_pdb = protein_input, 
#     chain_a = chain_1, 
#     chain_b = chain_2, 
#     return_breakdown = True
# )
# results_dict = {"Entropy type": "Free Energy contribution (J/(mol))",
#     "Translational-Rotational": results.dS_trans_rot, 
#     "Hydrophobic": results.dS_hydrophobic, 
#     "Side-chain": results.dS_sidechain, 
#     "Backbone": results.dS_backbone, 
#     "Total": results.dS_total
# }
# st.table(results_dict)