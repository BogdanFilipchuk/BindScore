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
interaction_color = {                                   # Define colors for interaction types
        "hydrogen_bond": '#e6194b',
        "pi_pi_stacking": "#3cb44b",
        "disulfide_bond": "#bfef45",
        "metal_coordination": "#f58231",
        "salt_bridge": "#911eb4",
        "hydrophobic_contact": "#4363d8",
        "halogen_bond": "#f032e6",
        "dipole-dipole": "#42d4f4",
        "unclassified": "#999999",
}

def get_interaction_color(energy, norm):                       
    cmap = mcolors.LinearSegmentedColormap.from_list("bgr",  # Create a custom colormap for the binding energy differences
        [(0.0, "#4363d8"), (0.5, "#999999"), (1.0, "#e6194b")],
        N=512,
    )
    return mcolors.to_hex(cmap(norm(energy)))           # Map the energy to a color using the colormap and normalization

def get_binding_energy_summary(interactions: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(interactions)
    summary = (
        df.groupby("interaction type")["binding energy"]
        .sum()
        .reset_index()
        .rename(columns={"interaction type": "Interaction Type", "binding energy": "Total Interaction Enthalpy [kJ/mol]"})
        .sort_values("Total Interaction Enthalpy [kJ/mol]"))
    def highlight_row(row):
        color = interaction_color.get(row["Interaction Type"], "white")
        return [f"background-color: {color}"] * len(row)
    return summary.style.apply(highlight_row, axis=1) 

#####TITLE#####
st.title("BindScore.\nBinding Enthalpy Predictor")

#####TEMPERATURE SLIDER#####
#st.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### protein input and visualization
protein_input_name = st.text_input("Enter protein PDB RCSB database ID (e.g. 6PYH)", "6PYH")
protein_input = fetch_pdb(protein_input_name)           # input the protein PDB ID
protein_input_parsed = Protein_Structure(protein_input) # turn data into a Protein_Structure object
chains = protein_input_parsed.chains()                  # get the chains for interaction analysis
# residues_importance = calculate_residue_importance(protein_input)

### selectbox with the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["cartoon", "stick"])  # the selectbox

### visualization button and viewer
if st.button("Visualize"):
    with st.spinner(text="Loading...", width="content"):  # Show spinner while processing
        
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
col1,col2 =st.columns(2)                                # Create two columns for chain selection
with col1:
    chain_1 = st.selectbox("Select chain A for interaction analysis", chains)
with col2:
    chain_2 = st.selectbox("Select chain B for interaction analysis", chains)
coloring_mode = st.selectbox("Select interaction coloring mode", ["By interaction type", "By energy contribution (blue - most negative)"]) # Select coloring mode for interactions

### interaction calculation and gradient coloring  
if st.button("Analyze the interactions"):               # Button to start the interaction analysis
    with st.spinner(text="Calculating... This could take up to a few minutes.", width="content"):  # Show spinner while processing
        
        interactions_full = get_dataset_interaction_list(protein_input_name, chain_1, chain_2) # Calculate interactions and binding energy differences
        energy_values = []
        classified_residues = set()
        for interaction in interactions_full:
            if interaction["interaction type"] != "unclassified":
                classified_residues.add(interaction["atom1"]["residue_seq"])
                classified_residues.add(interaction["atom2"]["residue_seq"])
            energy_values.append(interaction["binding energy"])
        interactions_filtered = [                                               # filter the interactions to exclude unclassified from meaningful interactions
        i for i in interactions_full
        if i["interaction type"] != "unclassified"
        or (i["atom1"]["residue_seq"] not in classified_residues
            and i["atom2"]["residue_seq"] not in classified_residues)
        ]
        if min(energy_values) < 0 and max(energy_values) > 0:   # Two-slope only for both pos. and neg. values.
            norm = mcolors.TwoSlopeNorm(                        # Normalize the binding energy around zero
                vmin=min(energy_values), vcenter=0.0, vmax=max(energy_values))
        else:                                                   # If all values are positive, use a regular normalization
            norm = mcolors.Normalize(vmin=min(energy_values), vmax=max(energy_values))

        viewer = py3Dmol.view(width=750, height=750)            # Create viewer
        viewer.addModel(protein_input, "pdb")                   # PDB format

        viewer.setStyle({'cartoon': {'color': "white"}})                     # Reset all to white cartoon first to make the interaction colors stand out
        if coloring_mode == "By energy contribution (blue - most negative)":              # Color interactions by energy contribution using the custom colormap    
            for interaction in interactions_filtered:                                                          
                    viewer.setStyle({'resi': interaction["atom1"]["residue_seq"]}, {'cartoon': {'color': get_interaction_color(interaction["binding energy"], norm) }})
                    viewer.setStyle({'resi': interaction["atom2"]["residue_seq"]}, {'cartoon': {'color': get_interaction_color(interaction["binding energy"], norm) }})
        else:                                                                               # Color interactions by type
            for interaction in interactions_filtered:
                viewer.setStyle({'resi': interaction["atom1"]["residue_seq"]}, {'cartoon': {'color': interaction_color.get(interaction["interaction type"]) }})
                viewer.setStyle({'resi': interaction["atom2"]["residue_seq"]}, {'cartoon': {'color': interaction_color.get(interaction["interaction type"]) }})

        viewer.zoomTo()                                     # Zoom to fit the molecule
        html = viewer._make_html()                          # Render in Streamlit
        components.html(html, width=750, height=750)        # Size of the streamlit component

### Interactions table
        enthalpy_data = get_binding_energy_summary(interactions_full)  # Prepare the interaction data for display in a table format
        st.text("The following interactions were identified.\nThe color of the interaction in the table corresponds to its color in the 3D view.")              # Title for the interactions table
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

# tests 

# interactions_full = get_dataset_interaction_list("6PYH", "A", "D")
# print([list([i["interaction type"], i["binding energy"], i["atom1"]["residue_seq"], i["atom2"]["residue_seq"]])   for i in interactions_sorted])