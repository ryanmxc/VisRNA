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
import pandas as pd
import scanpy as sc
from anndata import AnnData
from matplotlib.pyplot import rc_context

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

image = Image.open('./images/nsclc-logo.jpeg')

st.image(image, use_column_width=True)

st.write("""
# Cell type annotation and visualization.

""")
st.header('Cell type annotation')

#Load data
adata = st.session_state['adata']
import celltypist
from celltypist import models

# Indeed, the `model` argument defaults to `Immune_All_Low.pkl`.
model = models.Model.load(model = 'Human_Lung_Atlas.pkl')
predictions = celltypist.annotate(adata, model = 'Human_Lung_Atlas.pkl', majority_voting = True)
adata = predictions.to_adata()

sc.tl.umap(adata)

with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(adata, color = ['majority_voting'], legend_loc = 'on data')
    st.pyplot()