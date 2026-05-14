import json
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from pdbtime import fetch_pdb
from stmol import showmol
import py3Dmol

# def visualize_3D(molstring):
#     "Visualize the molecule in 3D using stmol"
#     w, h = 400, 400
#     xyzview = py3Dmol.view(width=w,height=w)
#     xyzview.addModel(molstring,'mol')
#     xyzview.setStyle({'sphere':{'colorscheme':'cyanCarbon', 'scale':0.25}, 'stick':{'colorscheme':'cyanCarbon'}})
#     xyzview.zoomTo()
#     xyzview.spin()
#     xyzview.setBackgroundColor('white')
#     showmol(xyzview, height = w,width=w)

#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####SIDEBAR#####
st.sidebar.slider("Temperature [K]", 0.0, 600.0, 0.7)

#####MAIN AREA#####
### choosing the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["skeleton", "cartoon", "surface"]) 

###protein input and visualization
protein_input = fetch_pdb(st.text_input("Enter protein PDB ID (e.g. 1A2B)", "1A2B"))  # input the protein PDB ID

if st.button("Load Protein"):
    html = f"""
    <script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

    <div id="viewer" style="width:100%; height:600px; position: relative;"></div>

    <script>
    const pdbData = {protein_input};

    const element = document.getElementById("viewer");

    const viewer = $3Dmol.createViewer(element, {{
        backgroundColor: "black"
    }});

    viewer.addModel(pdbData, "pdb");

    viewer.setStyle({{}}, {{
        cartoon: {{
            color: "spectrum"
        }}
    }});

    viewer.zoomTo();
    viewer.render();
    </script>
    """





    # # Create viewer
    # viewer = py3Dmol.view(width=900, height=700)
    # viewer.addModel(protein_input, "pdb")

    # # Style options
    # viewer.setStyle({"cartoon": {"color": "spectrum"}})

    # # Optional extras
    # viewer.addSurface(py3Dmol.VDW, {"opacity": 0.6})
    # viewer.zoomTo()

    # # Render in Streamlit
    # html = viewer._make_html()

    st.html(html, width=900, unsafe_allow_javascript=True)

### pH Slider (expander)
with st.expander("pH Scale"):
    fig, ax = plt.subplots(figsize=(10, 0.2)) # Create a wide, short figure (horizontal line) for the pH scale
    gradient = np.linspace(0, 1, 256).reshape(1, -1) # Create gradient
    ax.imshow(  #Show the gradient in graph
        gradient,
        aspect='auto',
        cmap='RdYlBu',
        extent=[0.0, 14.0, 0.0, 1.0]
    )
    ph = st.slider("Choose pH", 0.0, 14.0, 7.0)
    # Labels
    ax.set_yticks([])
    ax.set_xticks(range(15))
    ax.axvline(ph, color='black', linewidth=2) # Indicator line
    st.pyplot(fig)
    # pH Interpretation
    if ph < 6:
        st.error("Acidic")
    elif 6 <= ph < 7:
        st.success("Neutral")
    else:
        st.info("Basic")