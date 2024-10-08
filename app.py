import streamlit as st
import pandas as pd
import altair as alt

# Define molecules and their properties
molecules = [
    {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "properties": {
        "lipophilicity": 1.19, "solubility": 3, "permeability": 2, "clearance": 1, "potency": 2
    }},
    {"name": "Ibuprofen", "smiles": "CC(C)CC1=CC=C(C=C1)[C@H](C)C(=O)O", "properties": {
        "lipophilicity": 3.97, "solubility": 2, "permeability": 3, "clearance": 2, "potency": 2
    }},
    {"name": "Paracetamol", "smiles": "CC(=O)NC1=CC=C(C=C1)O", "properties": {
        "lipophilicity": 0.46, "solubility": 1, "permeability": 1, "clearance": 2, "potency": 1
    }},
    {"name": "Olaparib", "smiles": "C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=C(C=C3)CN4C(=O)C=CC4=O)F", "properties": {
        "lipophilicity": 1.95, "solubility": 2, "permeability": 3, "clearance": 2, "potency": 3
    }},
    {"name": "Abiraterone", "smiles": "CC12CCC3C(C1CCC2(C#C)O)CCC4C3(CCC(C4)O)C", "properties": {
        "lipophilicity": 5.19, "solubility": 1, "permeability": 3, "clearance": 1, "potency": 3
    }}
]

property_descriptions = {
    "lipophilicity": "Measure of a drug's ability to dissolve in fats, oils, and non-polar solvents",
    "solubility": "Ability of a drug to dissolve in water",
    "permeability": "Ease with which a drug can pass through cell membranes",
    "clearance": "Rate at which a drug is removed from the body",
    "potency": "Amount of drug required to produce a specific effect"
}

def get_traffic_light_color(value):
    if value <= 1:
        return "red"
    elif value <= 2:
        return "yellow"
    else:
        return "green"

st.title("Drug Property Viewer")

# Molecule selection
selected_molecules = st.multiselect("Select molecules (up to 5)", [m["name"] for m in molecules], default=["Aspirin"], max_selections=5)

# View toggle
view_type = st.radio("Select view type", ["Traffic Light", "Radar Plot"])

# Display properties
if view_type == "Traffic Light":
    for prop in property_descriptions.keys():
        st.subheader(prop)
        cols = st.columns(len(selected_molecules))
        for i, mol_name in enumerate(selected_molecules):
            mol = next(m for m in molecules if m["name"] == mol_name)
            value = mol["properties"][prop]
            color = get_traffic_light_color(value)
            cols[i].markdown(f"<h3 style='text-align: center; color: {color};'>{mol_name}</h3>", unsafe_allow_html=True)
            cols[i].markdown(f"<div style='width: 50px; height: 50px; border-radius: 25px; background-color: {color}; margin: auto;'></div>", unsafe_allow_html=True)

else:  # Radar Plot
    df = pd.DataFrame([m["properties"] for m in molecules if m["name"] in selected_molecules])
    df.index = selected_molecules

    # Prepare data for Altair
    df_melted = df.reset_index().melt(id_vars='index', var_name='property', value_name='value')
    df_melted = df_melted.rename(columns={'index': 'molecule'})

    # Create Radar Plot using Altair
    base = alt.Chart(df_melted).encode(
        theta=alt.Theta('property:N', sort=None),
        radius=alt.Radius('value:Q', scale=alt.Scale(type='sqrt', zero=True, rangeMin=20)),
        color='molecule:N'
    )
    
    chart = base.mark_line(closed=True).encode(
        alt.OpacityValue(0.2)
    ) + base.mark_point().encode(
        alt.OpacityValue(0.8)
    )
    
    st.altair_chart(chart, use_container_width=True)

# Property descriptions
st.subheader("Property Descriptions")
for prop, desc in property_descriptions.items():
    st.markdown(f"**{prop}**: {desc}")

# Display molecular structures
st.subheader("Molecular Structures")
for mol_name in selected_molecules:
    mol = next(m for m in molecules if m["name"] == mol_name)
    st.image(f"https://cactus.nci.nih.gov/chemical/structure/{mol['smiles']}/image", caption=mol_name)
