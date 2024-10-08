import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

st.set_page_config(page_title="Molecular Property Viewer", layout="wide")

@st.cache_data
def load_molecules():
    # This function will now cache the molecules data
    return [
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
    # Add 7 more molecules here to reach a total of 12
]

@st.cache_data
def mol_to_svg(smiles, size=150):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        rdDepictor.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(size, size)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')
    else:
        return "Invalid SMILES"

@st.cache_data
def get_property_descriptions():
    return {
        "lipophilicity": "Measure of a drug's ability to dissolve in fats, oils, and non-polar solvents",
        "solubility": "Ability of a drug to dissolve in water",
        "permeability": "Ease with which a drug can pass through cell membranes",
        "clearance": "Rate at which a drug is removed from the body",
        "potency": "Amount of drug required to produce a specific effect"
    }

def molecule_selection_page():
    st.subheader("Select Molecules (up to 4)")
    molecules = load_molecules()

    for i in range(0, len(molecules), 4):
        cols = st.columns(4)
        for j, col in enumerate(cols):
            if i + j < len(molecules):
                molecule = molecules[i + j]
                with col:
                    selected = st.checkbox("", key=f"mol_{i+j}", value=molecule['name'] in st.session_state.selected_molecules)
                    svg = mol_to_svg(molecule['smiles'])
                    if svg != "Invalid SMILES":
                        st.components.v1.html(svg, height=150, width=150)
                    else:
                        st.warning(f"Could not render molecule: {molecule['name']}")
                    st.markdown(f"<p style='text-align: center; font-weight: bold; font-size: 18px;'>{molecule['name']}</p>", unsafe_allow_html=True)
                
                if selected and molecule['name'] not in st.session_state.selected_molecules:
                    if len(st.session_state.selected_molecules) < 4:
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

    molecules = load_molecules()
    selected_data = [m for m in molecules if m['name'] in st.session_state.selected_molecules]

    # Display selected molecules
    st.subheader("Selected Molecules")
    cols = st.columns(4)
    for i, mol in enumerate(selected_data):
        with cols[i]:
            svg = mol_to_svg(mol['smiles'])
            if svg != "Invalid SMILES":
                st.components.v1.html(svg, height=150, width=150)
            else:
                st.warning(f"Could not render molecule: {mol['name']}")
            st.markdown(f"<p style='text-align: center; font-weight: bold; font-size: 18px;'>{mol['name']}</p>", unsafe_allow_html=True)

    view_type = st.radio("Select view type", ["Traffic Light", "Radar Plot"])

    if view_type == "Traffic Light":
        display_traffic_light(selected_data)
    else:
        display_radar_plot(selected_data)

    display_property_descriptions()

def display_traffic_light(selected_data):
    df = pd.DataFrame([{k: f"{v:.2f}" for k, v in m['properties'].items()} for m in selected_data])
    df.index = [m['name'] for m in selected_data]
    
    def color_cells(val, prop):
        color = get_traffic_light_color(prop, float(val))
        return f'background-color: {color}; color: black; text-align: center; vertical-align: middle;'

    styled_df = df.style.apply(lambda col: [color_cells(val, col.name) for val in col], axis=0)
    
    styled_df = styled_df.set_table_styles([
        {'selector': 'th', 'props': [('font-weight', 'bold'), ('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': 'td', 'props': [('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': '.col_heading', 'props': [('text-align', 'center'), ('font-weight', 'bold')]},
        {'selector': '.row_heading', 'props': [('text-align', 'center'), ('font-weight', 'bold')]},
    ])
    
    # Convert to HTML and adjust cell padding
    html = styled_df.to_html()
    html = html.replace('<td', '<td style="padding: 10px;"')
    html = html.replace('<th', '<th style="padding: 10px;"')
    
    # Wrap the table in a div with custom CSS for center alignment
    centered_html = f"""
    <div style="display: flex; justify-content: center;">
        <style>
            table {{
                border-collapse: collapse;
                margin: 0 auto;
            }}
            th, td {{
                border: 1px solid black;
                text-align: center;
                vertical-align: middle;
            }}
        </style>
        # {html}
    </div>
    """
    
    st.markdown(centered_html, unsafe_allow_html=True)
    
@st.cache_data
def prepare_radar_data(selected_data):
    df = pd.DataFrame([m["properties"] for m in selected_data])
    df = df.round(2)  # Round to 2 decimal places
    df.index = [m["name"] for m in selected_data]
    return df
    
def display_radar_plot(selected_data):
    df = prepare_radar_data(selected_data)
    
    # Define a color palette with distinct colors
    color_palette = px.colors.qualitative.Bold

    fig = go.Figure()

    for i, molecule in enumerate(df.index):
        color = color_palette[i % len(color_palette)]
        fig.add_trace(go.Scatterpolar(
            r=df.loc[molecule].values.tolist() + [df.loc[molecule].values[0]],
            theta=df.columns.tolist() + [df.columns[0]],
            fill='toself',
            name=molecule,
            line_color=color,
            fillcolor=color,
            opacity=0.5  # Adjust this value to change transparency (0.0 to 1.0)
        ))

    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 5],
                tickfont=dict(size=14)
            ),
            angularaxis=dict(
                tickfont=dict(size=14)
            )
        ),
        showlegend=True,
        legend=dict(font=dict(size=16)),
        title=dict(
            text="Molecular Properties Comparison",
            font=dict(size=20)
        ),
        height=600,
        width=800
    )

    st.plotly_chart(fig, use_container_width=True)

def display_property_descriptions():
    st.subheader("Property Descriptions")
    st.markdown(
        "<style>div[data-testid='stMarkdownContainer'] ul { font-size: 14px; }</style>",
        unsafe_allow_html=True
    )
    property_descriptions = get_property_descriptions()  # Use the cached function
    for prop, desc in property_descriptions.items():
        st.markdown(f"- **{prop}**: {desc}")

def get_traffic_light_color(property_name, value):
    thresholds = {
        "lipophilicity": {"low": 2, "high": 4},
        "solubility": {"low": 1, "high": 2},
        "permeability": {"low": 1, "high": 2},
        "clearance": {"low": 1, "high": 2},
        "potency": {"low": 1, "high": 2}
    }
    
    low, high = thresholds[property_name]["low"], thresholds[property_name]["high"]
    if value <= low:
        return "red"
    elif value <= high:
        return "yellow"
    else:
        return "green"

def main():
    st.title("Molecular Property Viewer")

    # Use session state to store selected molecules
    if 'selected_molecules' not in st.session_state:
        st.session_state.selected_molecules = []

    if 'page' not in st.session_state:
        st.session_state.page = 'selection'

    if st.session_state.page == 'selection':
        molecule_selection_page()
    elif st.session_state.page == 'property_view':
        property_view_page()

if __name__ == "__main__":
    main()
