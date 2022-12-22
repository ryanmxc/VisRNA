######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
import scanpy as sc
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import os
######################
# Page Title
######################

print(os.getcwd())
image = Image.open('./images/nsclc-logo.jpeg')

st.image(image, use_column_width=True)

st.write("""
## load your data one the left panel or use example file.
""")


######################
# Input scRNA-seq data (Side Panel)
######################

st.sidebar.header('User scRNA-seq data (please upload files)')
uploaded_file = st.sidebar.file_uploader("Choose a CSV file with scRNA-seq data", type="tsv")
use_example_file = st.sidebar.checkbox(
    "Use example file", False, key='1', help="Use in-built example file to demo the app"
)

# If CSV is not uploaded and checkbox is filled, use values from the example file
# and pass them down to the next if block
if use_example_file:
    uploaded_file = './Data/LUAD-003-01-1A.gene.expression.matrix.tsv'

st.header('Read scRNA-seq data using scanpy')

#Load data
if uploaded_file is not None:
    print(uploaded_file)
    st.write("Using uploaded file")
    ge_matrix = pd.read_csv(uploaded_file,sep='\t')
    st.write(ge_matrix)
    # ge_matrix = uploaded_file
else:
    st.write("Example file was preloaded")
    ge_matrix = './Data/LUAD-003-01-1A.gene.expression.matrix.tsv'
adata = sc.read(ge_matrix, cache=True)
adata = adata.transpose()
cell_count = adata.X.shape[0]
gene_count = adata.X.shape[1]
ge_df = pd.read_csv(ge_matrix,sep='\t')
ge_df.head()
st.write('Number of cells in the data: ', cell_count)
st.write('Number of genes in the data: ', gene_count)
st.session_state['adata'] = adata
st.session_state['ge_df'] = ge_df