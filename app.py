import streamlit as st
import pandas as pd
import altair as alt
from math import pi

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
    }},
    # Add 10 more molecules here to reach a total of 15
]

property_descriptions = {
    "lipophilicity": "Measure of a drug's ability to dissolve in fats, oils, and non-polar solvents",
    "solubility": "Ability of a drug to dissolve in water",
    "permeability": "Ease with which a drug can pass through cell membranes",
    "clearance": "Rate at which a drug is removed from the body",
    "potency": "Amount of drug required to produce a specific effect"
}

# Thresholds for traffic light coloring (easily adjustable)
thresholds = {
    "lipophilicity": {"low": 2, "high": 4},
    "solubility": {"low": 1, "high": 2},
    "permeability": {"low": 1, "high": 2},
    "clearance": {"low": 1, "high": 2},
    "potency": {"low": 1, "high": 2}
}

def get_traffic_light_color(property_name, value):
    low, high = thresholds[property_name]["low"], thresholds[property_name]["high"]
    if value <= low:
        return "red"
    elif value <= high:
        return "yellow"
    else:
        return "green"

st.set_page_config(page_title="Molecular Property Viewer", layout="wide")

# Initialize session state
if 'page' not in st.session_state:
    st.session_state.page = 'selection'
if 'selected_molecules' not in st.session_state:
    st.session_state.selected_molecules = []

def main():
    st.title("Molecular Property Viewer")

    if st.session_state.page == 'selection':
        molecule_selection_page()
    elif st.session_state.page == 'property_view':
        property_view_page()

def molecule_selection_page():
    st.subheader("Select Molecules (up to 5)")

    # Determine number of columns based on orientation
    cols = st.columns(5 if st.session_state.get('is_landscape', True) else 3)
    
    for i, molecule in enumerate(molecules):
        with cols[i % len(cols)]:
            selected = st.checkbox(molecule['name'], key=f"mol_{i}")
            st.image(f"https://cactus.nci.nih.gov/chemical/structure/{molecule['smiles']}/image", width=150)
            if selected and molecule['name'] not in st.session_state.selected_molecules:
                if len(st.session_state.selected_molecules) < 5:
                    st.session_state.selected_molecules.append(molecule['name'])
            elif not selected and molecule['name'] in st.session_state.selected_molecules:
                st.session_state.selected_molecules.remove(molecule['name'])

    if st.button("View Properties") and st.session_state.selected_molecules:
        st.session_state.page = 'property_view'

def property_view_page():
    if st.button("← Back to Selection"):
        st.session_state.page = 'selection'
        st.rerun()

    view_type = st.radio("Select view type", ["Traffic Light", "Radar Plot"])

    selected_data = [m for m in molecules if m['name'] in st.session_state.selected_molecules]

    if view_type == "Traffic Light":
        display_traffic_light(selected_data)
    else:
        display_radar_plot(selected_data)

    # Display molecular structures
    st.subheader("Molecular Structures")
    for mol in selected_data:
        st.image(f"https://cactus.nci.nih.gov/chemical/structure/{mol['smiles']}/image", caption=mol['name'])

def display_traffic_light(selected_data):
    for prop in property_descriptions.keys():
        st.subheader(prop)
        cols = st.columns(len(selected_data))
        for i, mol in enumerate(selected_data):
            value = mol["properties"][prop]
            color = get_traffic_light_color(prop, value)
            font_size = max(10, int(20 / len(selected_data)))  # Adjust font size based on number of molecules
            cols[i].markdown(f"<h3 style='text-align: center; color: {color}; font-size: {font_size}px;'>{mol['name']}</h3>", unsafe_allow_html=True)
            cols[i].markdown(f"<div style='width: 50px; height: 50px; border-radius: 25px; background-color: {color}; margin: auto;'></div>", unsafe_allow_html=True)

def display_radar_plot(selected_data):
    df = pd.DataFrame([m["properties"] for m in selected_data])
    df.index = [m["name"] for m in selected_data]
    df_melted = df.reset_index().melt(id_vars='index', var_name='property', value_name='value')
    df_melted = df_melted.rename(columns={'index': 'molecule'})

    chart = alt.Chart(df_melted).transform_calculate(
        angle=f"datum.property === '{df_melted['property'].iloc[-1]}' ? 0 : 2 * PI * (datum.property_index + 1) / {len(df_melted['property'].unique())}"
    ).mark_line(point=True).encode(
        x=alt.X('x:Q', axis=None),
        y=alt.Y('y:Q', axis=None),
        color='molecule:N',
        order='property_index:Q',
        detail='molecule:N'
    ).transform_calculate(
        property_index="indexof(datum.property, datum.property)",
        x=f"datum.value * cos(2 * PI * datum.property_index / {len(df_melted['property'].unique())})",
        y=f"datum.value * sin(2 * PI * datum.property_index / {len(df_melted['property'].unique())})"
    )

    st.altair_chart(chart, use_container_width=True)

if __name__ == "__main__":
    main()
