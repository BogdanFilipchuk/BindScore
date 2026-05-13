import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

st.title("Binding Energy Predictor")

#choosing the visualization type of the protein
Protein_drawmode = st.selectbox("Protein visualization mode", ["Skeleton", "Blah", "Blah"]) 

st.sidebar.slider("Temperature [K]", 0.0, 600.0, 0.7)

st.write("Main content area")

# pH Slider
ph = st.slider("pH", 0.0, 14.0, 0.5)

fig, ax = plt.subplots(figsize=(10, 1)) # Create gradient
gradient = np.linspace(0, 1, 256).reshape(1, -1)
ax.imshow(
    gradient,
    aspect='auto',
    cmap='RdYlGn',
    extent=[0, 14, 0, 1]
)

ax.axvline(ph, color='black', linewidth=3) # Indicator line

# Labels
ax.set_yticks([])
ax.set_xticks(range(15))
ax.set_xlabel("pH Scale")

st.pyplot(fig)

# Description
if ph < 6.5:
    st.error("Acidic")
elif 6.5 <= ph < 7.5:
    st.success("Neutral")
else:
    st.info("Basic / Alkaline")