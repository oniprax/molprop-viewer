import base64
import io
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

st.set_page_config(page_title="Molecular Property Predictor", layout="wide")

def local_css(file_name):
    with open(file_name, 'r') as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

def display_large_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(400, 400))
    st.image(img, use_column_width=True)
    
@st.cache_data
def load_molecule_dataframe():
    df = pd.read_pickle("./ccdd_smidf.pkl")
    # Convert SMILES to RDKit molecule objects
    df['Molecule'] = df['Mol'].apply(Chem.MolFromSmiles)
    df['R1'] = df['R1'].apply(Chem.MolFromSmiles)
    df['R2'] = df['R2'].apply(Chem.MolFromSmiles)
    df['R3'] = df['R3'].apply(Chem.MolFromSmiles)
    df['R4'] = df['R4'].apply(Chem.MolFromSmiles)
    return df


def mol_to_img(mol):
    img = Draw.MolToImage(mol, size=(150, 150))
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def display_molecule_table(df):
    # Initialize session state for selections if not exists
    if 'selections' not in st.session_state:
        st.session_state.selections = [False] * len(df)

    # Create a column for checkboxes
    df['Select'] = st.session_state.selections

    # Function to create HTML for molecule image
    def mol_to_html(mol):
        return f'<img src="data:image/png;base64,{mol_to_img(mol)}" width="150">'

    # Apply the function to create a new column with HTML
    df['Structure'] = df['Molecule'].apply(mol_to_html)

    # Reorder columns
    columns = ['Select', 'ID', 'Structure', 'R1', 'R2', 'R3', 'R4']
    df = df[columns]

    # Display the dataframe
    edited_df = st.data_editor(
        df,
        hide_index=True,
        column_config={
            "Select": st.column_config.CheckboxColumn(required=True),
            "Structure": st.column_config.Column(width="medium"),
            "ID": st.column_config.TextColumn(width="small"),
            "R1": st.column_config.TextColumn(width="medium"),
            "R2": st.column_config.TextColumn(width="medium"),
            "R3": st.column_config.TextColumn(width="medium"),
            "R4": st.column_config.TextColumn(width="medium"),
        },
        disabled=df.columns.drop('Select'),
        key="molecule_table"
    )

    # Update session state with new selections
    st.session_state.selections = edited_df['Select'].tolist()

    # Get selected molecules
    selected_molecules = edited_df[edited_df['Select']]['ID'].tolist()
    
    return selected_molecules
    
@st.cache_data
def load_molecules():
    # This function will now cache the molecules data
    return [{'name': 'ID01',
  'smiles': 'Cn1cc(cn1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.39,
                 'lipophilicity': 3.59,
                 'solubility': 39.0,
                 'permeability': 8.31,
                 'synthesisability': 0.72}},
 {'name': 'ID02',
  'smiles': 'CN1C(=O)C=C(c2cnn(c2)C)c2cc(ccc21)Nc1ccncc1Cl',
  'properties': {'potency': 4.91,
                 'lipophilicity': 3.73,
                 'solubility': 0.0,
                 'permeability': 10.51,
                 'synthesisability': 0.74}},
 {'name': 'ID03',
  'smiles': 'Clc1nccc(c1C#N)Nc1ccc2c(c1)C(=CC(=O)N2)c1cnn(c1)CCN1CCOCC1',
  'properties': {'potency': 5.43,
                 'lipophilicity': 3.39,
                 'solubility': 52.0,
                 'permeability': 6.85,
                 'synthesisability': 0.71}},
 # {'name': 'ID',
 #  'smiles': 'CN1CCN(CC1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
 #  'properties': {'potency': 5.44,
 #                 'lipophilicity': 2.94,
 #                 'solubility': 40.0,
 #                 'permeability': 5.52,
 #                 'synthesisability': 0.73}},
 {'name': 'ID04',
  'smiles': 'Cn1cc(cn1)C1=CC(=O)N(CC2CCOCC2)c2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.73,
                 'lipophilicity': 4.49,
                 'solubility': 86.0,
                 'permeability': 15.7,
                 'synthesisability': 0.7}},
 {'name': 'ID05',
  'smiles': 'CNC(=O)Cn1cc(cn1)C1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.82,
                 'lipophilicity': 2.8,
                 'solubility': 34.0,
                 'permeability': 4.25,
                 'synthesisability': 0.71}},
 {'name': 'ID06',
  'smiles': 'CNC(=O)[C@@H](C)NC1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 6.6,
                 'lipophilicity': 2.74,
                 'solubility': 64.0,
                 'permeability': 4.16,
                 'synthesisability': 0.67}},
 {'name': 'ID07',
  'smiles': 'CNC(=O)[C@H](C)NC1=CC(=O)Nc2ccc(cc21)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 5.57,
                 'lipophilicity': 2.74,
                 'solubility': 64.0,
                 'permeability': 4.16,
                 'synthesisability': 0.67}},
 # {'name': 'ID',
 #  'smiles': 'CCNC1=CC(=O)N(C)c2ccc(cc21)Nc1cc(nc(c1C#N)Cl)C(=O)O',
 #  'properties': {'potency': 6.05,
 #                 'lipophilicity': 3.33,
 #                 'solubility': 237.0,
 #                 'permeability': 6.37,
 #                 'synthesisability': 0.72}},
 {'name': 'ID08',
  'smiles': 'C[C@@H](NC1=CC(=O)N(C)c2ccc(cc21)Nc1ccnc(c1C#N)Cl)C(=O)NC1CCCC1',
  'properties': {'potency': 6.59,
                 'lipophilicity': 4.06,
                 'solubility': 58.0,
                 'permeability': 10.97,
                 'synthesisability': 0.67}},
 {'name': 'ID09',
  'smiles': 'C[C@H]1COC(=O)C2=C(N1)c1cc(ccc1N(C)C2=O)Nc1ccnc(c1C#N)Cl',
  'properties': {'potency': 7.03,
                 'lipophilicity': 3.17,
                 'solubility': 0.0,
                 'permeability': 5.96,
                 'synthesisability': 0.64}},
 # {'name': 'ID',
 #  'smiles': 'CN1C(=O)C=C(NC(C)(C)c2ncccn2)c2cc(ccc21)Nc1ccnc(c1C#N)Cl',
 #  'properties': {'potency': 6.33,
 #                 'lipophilicity': 4.34,
 #                 'solubility': 15.0,
 #                 'permeability': 13.51,
 #                 'synthesisability': 0.69}},
 {'name': 'ID10',
  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)-n1nc(cc1C)C',
  'properties': {'potency': 7.95,
                 'lipophilicity': 4.75,
                 'solubility': 21.07,
                 'permeability': 18.67,
                 'synthesisability': 0.61}},
 {'name': 'ID11',
  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)N1CC(F)C1',
  'properties': {'potency': 8.21,
                 'lipophilicity': 4.1,
                 'solubility': 17.0,
                 'permeability': 12.6,
                 'synthesisability': 0.6}},
 # {'name': 'ID',
 #  'smiles': 'CN1C(=O)C2=C(N[C@@H](C3CC3)C(F)(F)CO2)c2cc(ccc21)Nc1nc(ncc1Cl)N1C2CCC1CN(C2)C(=O)C',
 #  'properties': {'potency': 9.04,
 #                 'lipophilicity': 4.14,
 #                 'solubility': 34.81,
 #                 'permeability': 11.96,
 #                 'synthesisability': 0.48}},
 {'name': 'ID12',
  'smiles': 'C[C@@H]1C[C@H](O)CN(C1)c1ncc(c(n1)Nc1ccc2c(c1)C1=C(OCC(F)(F)[C@@H](N1)C1CC1)C(=O)N2C)Cl',
  'properties': {'potency': 8.56,
                 'lipophilicity': 4.15,
                 'solubility': 7.69,
                 'permeability': 12.03,
                 'synthesisability': 0.55}
 }
           ]


@st.cache_data
def mol_to_svg(smiles, size=150):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        rdDepictor.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DSVG(size, size)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().additionalAtomLabelPadding = 0.1
        drawer.drawOptions().bondLineWidth = 2
        drawer.drawOptions().minFontSize = 8
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg.replace('svg:', '')
    else:
        return "Invalid SMILES"

def get_traffic_light_color(property_name, value):
    thresholds = {
        "potency": {"low": 6, "high": 7},
        "lipophilicity": {"low": 3, "high": 4},
        "solubility": {"low": 30, "high": 60},
        "permeability": {"low": 8, "high": 12},
        "synthesisability": {"low": 0.65, "high": 0.7}
    }
    
    low, high = thresholds[property_name]["low"], thresholds[property_name]["high"]
    if value <= low:
        return "red"
    elif value <= high:
        return "yellow"
    else:
        return "green"
        
@st.cache_data
def get_property_descriptions():
    return {
        "potency": "Amount of molecule required to produce a specific effect",
        "lipophilicity": "Measure of a molecule's ability to dissolve in fats, oils, and non-polar solvents",
        "solubility": "Ability of a molecule to dissolve in water",
        "permeability": "Ease with which a molecule can pass through cell membranes",
        "synthesisability": "Ease with which a molecule can be synthesised in the lab"
    }
    
def molecule_selection_page():
    st.title("Molecule Selection")
    
    # Display large molecule
    core_smiles = 'O=C1C([*:3])=C([*:2])C2=CC(N[*:1])=CC=C2N1[*:4]'
    display_large_molecule(core_smiles)
    
    # Load and display DataFrame
    df = load_molecule_dataframe()
    
    st.write("Select up to 4 molecules:")
    selected_molecules = display_molecule_table(df)
    
    if len(selected_molecules) > 4:
        st.warning("Please select no more than 4 molecules.")
    
    if st.button("View Properties", key='view_properties') and 0 < len(selected_molecules) <= 4:
        st.session_state.selected_molecules = selected_molecules
        st.session_state.page = 'property_view'
        st.rerun()

def property_view_page():
    if st.button("← Back to Selection"):
        st.session_state.page = 'selection'
        st.rerun()

    molecules = load_molecules()
    selected_data = [m for m in molecules if m['name'] in st.session_state.selected_molecules]

    # Custom CSS to reduce space between boxes and increase size
    st.markdown("""
        <style>
        .element-container {
            margin: 0px !important;
            padding: 0px !important;
        }
        </style>
    """, unsafe_allow_html=True)

    # Display selected molecules
    st.subheader("Selected Molecules")
    cols = st.columns(4)
    for i, mol in enumerate(selected_data):
        with cols[i]:
            svg = mol_to_svg(mol['smiles'], size=150)  # Increased size
            if svg != "Invalid SMILES":
                st.components.v1.html(svg, height=180, width=180)
            else:
                st.warning(f"Could not render molecule: {mol['name']}")
            st.markdown(f"<p style='text-align: center; font-weight: bold; font-size: 20px; margin: 0;'>{mol['name']}</p>", unsafe_allow_html=True)

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
    
    # Calculate the max value for each property
    max_values = df.max()
    overall_max = max(max_values) * 1.1  # Add 10% padding

    # Scale the data to fit within 0-1 range
    df_scaled = df / max_values

    fig = go.Figure()

    for i,molecule in enumerate(df_scaled.index):
        color = color_palette[i % len(color_palette)]
        fig.add_trace(go.Scatterpolar(
            r=df_scaled.loc[molecule].values.tolist() + [df_scaled.loc[molecule].values[0]],
            theta=df_scaled.columns.tolist() + [df_scaled.columns[0]],
            fill='toself',
            name=molecule,
            line_color=color,
            fillcolor=color,
            opacity=0.3
        ))

    # Update layout
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 1],  # Now all data is scaled to 0-1
                tickfont=dict(size=14),
                tickmode='array',
                tickvals=np.linspace(0, 1, 6),
                ticktext=[f'{v:.1f}' for v in np.linspace(0, overall_max, 6)]
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

    # Add property names as annotations
    # for i, prop in enumerate(df.columns):
    #     angle = i * 360 / len(df.columns)
    #     fig.add_annotation(
    #         # text=f"{prop}<br>Max: {max_values[prop]:.2f}",
    #         text=f"{prop}",
    #         x=1.2 * np.cos(np.radians(angle)),
    #         y=1.2 * np.sin(np.radians(angle)),
    #         showarrow=False,
    #         font=dict(size=12)
    #     )

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
        "potency": {"low": 6, "high": 7},
        "lipophilicity": {"low": 3, "high": 4},
        "solubility": {"low": 30, "high": 60},
        "permeability": {"low": 8, "high": 12},
        "synthesisability": {"low": 0.65, "high": 0.7}
    }
    
    low, high = thresholds[property_name]["low"], thresholds[property_name]["high"]
    if value <= low:
        return "red"
    elif value <= high:
        return "yellow"
    else:
        return "green"

def main():
    st.set_page_config(page_title="Molecular Property Predictor", layout="wide")

    if 'page' not in st.session_state:
        st.session_state.page = 'selection'

    if st.session_state.page == 'selection':
        molecule_selection_page()
    elif st.session_state.page == 'property_view':
        property_view_page()

if __name__ == "__main__":
    main()
