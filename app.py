import streamlit as st
import pandas as pd
import altair as alt
from rdkit import Chem
from rdkit.Chem import Draw
import io
import base64
import math

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
    # Add more molecules here
]

property_descriptions = {
    "lipophilicity": "Measure of a drug's ability to dissolve in fats, oils, and non-polar solvents",
    "solubility": "Ability of a drug to dissolve in water",
    "permeability": "Ease with which a drug can pass through cell membranes",
    "clearance": "Rate at which a drug is removed from the body",
    "potency": "Amount of drug required to produce a specific effect"
}

# Thresholds for traffic light coloring
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

def mol_to_img(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode()

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

    cols = st.columns(5)
    for i, molecule in enumerate(molecules):
        with cols[i % 5]:
            selected = st.checkbox(molecule['name'], key=f"mol_{i}")
            st.image(f"data:image/png;base64,{mol_to_img(molecule['smiles'])}", width=150)
            if selected and molecule['name'] not in st.session_state.selected_molecules:
                if len(st.session_state.selected_molecules) < 5:
                    st.session_state.selected_molecules.append(molecule['name'])
            elif not selected and molecule['name'] in st.session_state.selected_molecules:
                st.session_state.selected_molecules.remove(molecule['name'])

    if st.button("View Properties", key='view_properties') and st.session_state.selected_molecules:
        st.session_state.page = 'property_view'
        st.rerun()

def property_view_page():
    if st.button("← Back to Selection"):
        st.session_state.page = 'selection'
        st.rerun()

    selected_data = [m for m in molecules if m['name'] in st.session_state.selected_molecules]

    # Display selected molecules
    st.subheader("Selected Molecules")
    cols = st.columns(len(selected_data))
    for i, mol in enumerate(selected_data):
        with cols[i]:
            st.image(f"data:image/png;base64,{mol_to_img(mol['smiles'])}", caption=mol['name'], width=150)

    view_type = st.radio("Select view type", ["Traffic Light", "Radar Plot"])

    if view_type == "Traffic Light":
        display_traffic_light(selected_data)
    else:
        display_radar_plot(selected_data)

def display_traffic_light(selected_data):
    df = pd.DataFrame([m['properties'] for m in selected_data])
    df.index = [m['name'] for m in selected_data]
    
    # Create a styled dataframe
    def color_cells(val, prop):
        color = get_traffic_light_color(prop, val)
        return f'background-color: {color}; color: black; font-weight: bold; text-align: center; vertical-align: middle;'

    styled_df = df.style.apply(lambda col: [color_cells(val, col.name) for val in col], axis=0)
    
    # Center-align and bold the column names and index
    styled_df = styled_df.set_table_styles([
        {'selector': 'th', 'props': [('font-weight', 'bold'), ('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': 'td', 'props': [('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': 'tr:hover', 'props': [('background-color', 'lightgrey')]},
    ])
    
    # Bold and center the index (molecule names)
    styled_df = styled_df.set_properties(**{'font-weight': 'bold', 'text-align': 'center', 'vertical-align': 'middle'})
    
    # Convert to HTML and adjust cell padding
    html = styled_df.to_html()
    html = html.replace('<td', '<td style="padding: 10px;"')
    
    st.write(html, unsafe_allow_html=True)

def display_radar_plot(selected_data):
    df = pd.DataFrame([m["properties"] for m in selected_data])
    df.index = [m["name"] for m in selected_data]
    
    # Melt the dataframe
    df_melted = df.reset_index().melt(id_vars='index', var_name='property', value_name='value')
    df_melted = df_melted.rename(columns={'index': 'molecule'})

    # Calculate angles for each property
    properties = list(df.columns)
    n_properties = len(properties)
    angles = [i / n_properties * 2 * math.pi for i in range(n_properties)]
    angles += angles[:1]  # Repeat the first angle to close the polygon

    # Create the base chart
    base = alt.Chart(df_melted).encode(
        theta=alt.Theta('property:N', sort=None, stack=None),
        radius=alt.Radius('value:Q', scale=alt.Scale(type='sqrt', zero=True, rangeMin=20)),
        color='molecule:N'
    )

    # Create the radar chart
    lines = base.mark_line(point=True).encode(
        order='property'
    )

    # Add circular grid lines
    grid_data = pd.DataFrame({'radius': [1, 2, 3, 4, 5]})
    grid = alt.Chart(grid_data).encode(
        theta=alt.Theta(datum=2 * math.pi, scale=alt.Scale(domain=[0, 2 * math.pi])),
        radius=alt.Radius('radius:Q', scale=alt.Scale(domain=[0, 5]))
    ).mark_circle(
        color='lightgray',
        strokeWidth=1,
        opacity=0.5,
        fill=None
    )

    # Add property labels
    labels = alt.Chart(pd.DataFrame({
        'property': properties,
        'angle': angles[:-1],
        'radius': [5.5] * n_properties
    })).encode(
        text='property:N',
        theta='angle:Q',
        radius='radius:Q'
    ).mark_text(align='center', baseline='middle')

    # Combine all elements
    radar = (grid + lines + labels).properties(width=500, height=500)

    st.altair_chart(radar, use_container_width=True)

if __name__ == "__main__":
    main()
