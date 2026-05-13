import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

#####TITLE#####
st.title("BindScore - Binding Energy Predictor")

#####MAIN AREA#####

### choosing the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["Skeleton", "Blah", "Blah"]) 

st.sidebar.slider("Temperature [K]", 0.0, 600.0, 0.7)

### pH Slider
fig, ax = plt.subplots(figsize=(10, 0.2)) # Create a wide, short figure (horizontal line) for the pH scale
gradient = np.linspace(0, 1, 256).reshape(1, -1) # Create gradient
ax.imshow(  #Show the gradient in graph
    gradient,
    aspect='auto',
    cmap='RdYlBu',
    extent=[0.0, 14.0, 0.0, 1.0]
)
with st.expander("pH Scale"):
    ph = st.slider("Choose pH", 0.0, 14.0, 7.0)
    # Labels
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



