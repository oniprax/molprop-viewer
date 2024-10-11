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
    return [{'name': 'ID01',
  'smiles': 'Cn1cc(cn1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.39,
                 'lipophilicity': 3.59,
                 'solubility': 39.0,
                 'permeability': 8.31,
                 'synthetic accessibility': 0.72}},
 {'name': 'ID02',
  'smiles': 'CN1C(=O)C=C(c2cnn(c2)C)c2cc(ccc21)Nc1ccncc1Cl',
  'properties': {'potency': 4.91,
                 'lipophilicity': 3.73,
                 'solubility': 0.0,
                 'permeability': 10.51,
                 'synthetic accessibility': 0.74}},
 {'name': 'ID03',
  'smiles': 'Clc1nccc(c1C#N)Nc1ccc2c(c1)C(=CC(=O)N2)c1cnn(c1)CCN1CCOCC1',
  'properties': {'potency': 5.43,
                 'lipophilicity': 3.39,
                 'solubility': 52.0,
                 'permeability': 6.85,
                 'synthetic accessibility': 0.71}},
 # {'name': 'ID04',
 #  'smiles': 'CN1CCN(CC1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
 #  'properties': {'potency': 5.44,
 #                 'lipophilicity': 2.94,
 #                 'solubility': 40.0,
 #                 'permeability': 5.52,
 #                 'synthetic accessibility': 0.73}},
 {'name': 'ID05',
  'smiles': 'Cn1cc(cn1)C1=CC(=O)N(CC2CCOCC2)c2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.73,
                 'lipophilicity': 4.49,
                 'solubility': 86.0,
                 'permeability': 15.7,
                 'synthetic accessibility': 0.7}},
 {'name': 'ID06',
  'smiles': 'CNC(=O)Cn1cc(cn1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.82,
                 'lipophilicity': 2.8,
                 'solubility': 34.0,
                 'permeability': 4.25,
                 'synthetic accessibility': 0.71}},
 {'name': 'ID07',
  'smiles': 'CNC(=O)[C@@H](C)NC1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 6.6,
                 'lipophilicity': 2.74,
                 'solubility': 64.0,
                 'permeability': 4.16,
                 'synthetic accessibility': 0.67}},
 {'name': 'ID08',
  'smiles': 'CNC(=O)[C@H](C)NC1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.57,
                 'lipophilicity': 2.74,
                 'solubility': 64.0,
                 'permeability': 4.16,
                 'synthetic accessibility': 0.67}},
 # {'name': 'ID09',
 #  'smiles': 'CCNC1=CC(=O)N(C)c2ccc(cc21)Nc1cc(nc(c1C#N)Cl)C(=O)O',
 #  'properties': {'potency': 6.05,
 #                 'lipophilicity': 3.33,
 #                 'solubility': 237.0,
 #                 'permeability': 6.37,
 #                 'synthetic accessibility': 0.72}},
 {'name': 'ID10',
  'smiles': 'C[C@@H](NC1=CC(=O)N(C)c2ccc(cc21)Nc1ccnc(c1C#N)Cl)C(=O)NC1CCCC1',
  'properties': {'potency': 6.59,
                 'lipophilicity': 4.06,
                 'solubility': 58.0,
                 'permeability': 10.97,
                 'synthetic accessibility': 0.67}},
 {'name': 'ID11',
  'smiles': 'C[C@H]1COC(=O)C2=C(N1)c1cc(ccc1N(C)C2=O)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 7.03,
                 'lipophilicity': 3.17,
                 'solubility': 0.0,
                 'permeability': 5.96,
                 'synthetic accessibility': 0.64}},
 # {'name': 'ID12',
 #  'smiles': 'CN1C(=O)C=C(NC(C)(C)c2ncccn2)c2cc(ccc21)Nc1ccnc(c1C#N)Cl',
 #  'properties': {'potency': 6.33,
 #                 'lipophilicity': 4.34,
 #                 'solubility': 15.0,
 #                 'permeability': 13.51,
 #                 'synthetic accessibility': 0.69}},
 {'name': 'ID13',
  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)-n1nc(cc1C)C',
  'properties': {'potency': 7.95,
                 'lipophilicity': 4.75,
                 'solubility': 21.07,
                 'permeability': 18.67,
                 'synthetic accessibility': 0.61}},
 {'name': 'ID14',
  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)N1CC(F)C1',
  'properties': {'potency': 8.21,
                 'lipophilicity': 4.1,
                 'solubility': 17.0,
                 'permeability': 12.6,
                 'synthetic accessibility': 0.6}},
 # {'name': 'ID15',
 #  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)N1C2CCC1CN(C2)C(=O)C',
 #  'properties': {'potency': 9.04,
 #                 'lipophilicity': 4.14,
 #                 'solubility': 34.81,
 #                 'permeability': 11.96,
 #                 'synthetic accessibility': 0.48}},
 {'name': 'ID16',
  'smiles': 'C[C@@H]1C[C@H](O)CN(C1)c1ncc(c(n1)Nc1ccc2c(c1)C1=C(OCC(F)(F)[C@@H](N1)C1CC1)C(=O)N2C)Cl',
  'properties': {'potency': 8.56,
                 'lipophilicity': 4.15,
                 'solubility': 7.69,
                 'permeability': 12.03,
                 'synthetic accessibility': 0.55}
 }
           ]


@st.cache_data
def mol_to_svg(smiles, size=150):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        rdDepictor.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(size, size)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().additionalAtomLabelPadding = 0.3
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

    # Custom CSS to reduce padding
    st.markdown("""
        <style>
        .stSelectbox, .stCheckbox {
            padding-bottom: 0px;
        }
        .element-container {
            margin-bottom: 10px;
        }
        </style>
    """, unsafe_allow_html=True)

    # Use st.columns with custom widths
    col1, col2, col3, col4 = st.columns([1, 1, 1, 1])
    columns = [col1, col2, col3, col4]

    for i, molecule in enumerate(molecules):
        with columns[i % 4]:
            selected = st.checkbox("", key=f"mol_{i}", value=molecule['name'] in st.session_state.selected_molecules)
            svg = mol_to_svg(molecule['smiles'], size=150)  # Change size
            if svg != "Invalid SMILES":
                st.components.v1.html(svg, height=120, width=120)
            else:
                st.warning(f"Could not render molecule: {molecule['name']}")
            st.markdown(f"<p style='text-align: center; font-weight: bold; font-size: 14px; margin: 0;'>{molecule['name']}</p>", unsafe_allow_html=True)
        
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

    # Custom CSS to reduce padding
    st.markdown("""
        <style>
        .element-container {
            margin-bottom: 10px;
        }
        </style>
    """, unsafe_allow_html=True)

    # Display selected molecules
    st.subheader("Selected Molecules")
    cols = st.columns(4)
    for i, mol in enumerate(selected_data):
        with cols[i]:
            svg = mol_to_svg(mol['smiles'], size=150)  # Change size
            if svg != "Invalid SMILES":
                st.components.v1.html(svg, height=120, width=120)
            else:
                st.warning(f"Could not render molecule: {mol['name']}")
            st.markdown(f"<p style='text-align: center; font-weight: bold; font-size: 14px; margin: 0;'>{mol['name']}</p>", unsafe_allow_html=True)

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
        return f'background-color: {color}; color: black;'

    styled_df = df.style.apply(lambda col: [color_cells(val, col.name) for val in col], axis=0)
    
    styled_df = styled_df.set_table_styles([
        {'selector': 'th', 'props': [('font-weight', 'bold'), ('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': 'td', 'props': [('text-align', 'center'), ('vertical-align', 'middle')]},
        {'selector': '.col_heading', 'props': [('text-align', 'center'), ('font-weight', 'bold')]},
        {'selector': '.row_heading', 'props': [('text-align', 'center'), ('font-weight', 'bold')]},
    ])
    
    # Use Streamlit's native table function with custom CSS
    st.markdown("""
    <style>
        .stTable {
            width: 100%;
            text-align: center;
        }
        .stTable th {
            text-align: center !important;
        }
        .stTable td {
            text-align: center !important;
        }
    </style>
    """, unsafe_allow_html=True)
    
    st.table(styled_df)
    
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
