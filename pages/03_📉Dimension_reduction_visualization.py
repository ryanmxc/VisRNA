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

######################
# Page Title
######################

image = Image.open('./images/nsclc-logo.jpeg')

st.image(image, use_column_width=True)

st.write("""
# Dimension reduction and visualization.

""")
st.header('Dimension reduction visualization')

#Load data
adata = st.session_state['adata']
ge_df = st.session_state['ge_df']
gene_id = pd.DataFrame(columns=['gene'])
gene_id['gene'] = list(ge_df.index)
gene_id.head()
adata.var = gene_id

from matplotlib.pyplot import rc_context
st.set_option('deprecation.showPyplotGlobalUse', False)
st.subheader('Visualization by PCA')
sc.pp.pca(adata, n_comps=50, use_highly_variable=False, svd_solver='arpack')
with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.pca_scatter(adata)
    st.pyplot()

st.subheader('Visualization by UMAP')
# t-SNE
tsne_n_pcs = 20 # Number of principal components to use for t-SNE

# k-means
k = 18 # Number of clusters for k-means

# KNN
n_neighbors = 15 # Number of nearest neighbors for KNN graph
knn_n_pcs = 50 # Number of principal components to use for finding nearest neighbors

# UMAP
umap_min_dist = 0.3 
umap_spread = 1.0
# KNN graph
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs)
# UMAP
sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread)
# Louvain clustering
st.write('UMAP visualization with Louvain clustering')
sc.tl.louvain(adata)
with rc_context({'figure.figsize': (10, 10)}):
# Plot
    sc.pl.umap(adata, color=["louvain"],legend_loc='on data', legend_fontsize=10)
    st.pyplot()

# Leiden clustering
st.write('UMAP visualization with Leiden clustering')
sc.tl.leiden(adata)
# Plot
with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(adata, color=["leiden"],legend_loc='on data', legend_fontsize=10)
    st.pyplot()

st.subheader('Visualization by t-SNE')
sc.tl.tsne(adata, n_pcs=tsne_n_pcs)
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=18, random_state=0).fit(adata.obsm['X_pca'])
adata.obs['kmeans'] = kmeans.labels_.astype(str)
with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.tsne(adata, color=["kmeans"])
    st.pyplot()