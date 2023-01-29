######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
import scanpy as sc
import plotly.express as px

st.write("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Arimo');
html, body, [class*="css"]  {
   font-family: 'Arimo';
}
</style>
""", unsafe_allow_html=True)

######################
# Page Title
######################

image = Image.open('./images/load-data.png')
st.image(image, use_column_width=True)


######################
# Input scRNA-seq data (Side Panel)
######################

# st.sidebar.header('User scRNA-seq data (please upload files)')
uploaded_file = st.sidebar.file_uploader("Choose a CSV/TSV file with scRNA-seq data", type="tsv")
use_example_file = st.sidebar.checkbox(
    "Use example file", False, key='1', help="Use in-built example file to demo the app"
)


# If CSV is not uploaded and checkbox is filled, use values from the example file
# and pass them down to the next if block

st.header('Read scRNA-seq data using scanpy')

#Load data
if uploaded_file is not None:
    st.write("Using uploaded file")
    ge_df = pd.read_csv(uploaded_file,sep='\t')

elif use_example_file:
    st.write("An example file was loaded")
    ge_matrix = './Data/example_data.tsv'
    ge_df = pd.read_csv(ge_matrix,sep='\t',index_col=0)
 
else:
    st.write("An example file was preloaded")
    ge_matrix = './Data/example_data.tsv'
    ge_df = pd.read_csv(ge_matrix,sep='\t',index_col=0)

# Download example file
with open('./Data/example_data.tsv', 'rb') as f:
    s = f.read()
st.sidebar.download_button(
    label="Download example file",
    data=s,
    file_name='example_data.tsv',
    # mime='tsv',
)

st.header('Characteristics of loaded scRNA-seq data')

adata = sc.AnnData(ge_df)
adata = adata.transpose()
cell_count = adata.X.shape[0]
gene_count = adata.X.shape[1]
col1, col2 = st.columns(2)
col1.metric("Number of cells", cell_count)
col2.metric("Number of genes", gene_count)

st.session_state['adata'] = adata
st.session_state['ge_df'] = ge_df


example_genes = list(ge_df.index)[:10]
option = st.selectbox(
    'Select an example gene to show the distribution',
    example_genes)


st.write('You selected:', option)
if option:
    df = ge_df.T
    fig = px.histogram(df, x=option)
    st.plotly_chart(fig, use_container_width=True)