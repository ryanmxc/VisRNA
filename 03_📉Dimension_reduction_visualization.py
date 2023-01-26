######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
import scanpy as sc
import scipy as sp
from matplotlib import rcParams
import seaborn as sb
from sklearn.cluster import KMeans

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

image = Image.open('./images/dimension-reduction.png')

st.image(image, use_column_width=True)

st.subheader('Perform dimension reduction using scRNA-seq data.')

#Load data
adata = st.session_state['adata']
ge_df = st.session_state['ge_df']
gene_id = pd.DataFrame(columns=['gene'])
gene_id['gene'] = list(ge_df.index)
gene_id.head()
adata.var = gene_id
adata.var_names = gene_id['gene'].tolist()

from matplotlib.pyplot import rc_context
st.set_option('deprecation.showPyplotGlobalUse', False)
st.subheader('Visualization by PCA')

st.markdown('<div style="text-align: justify;"> PCA is a dimension reduction method \
</div>', unsafe_allow_html=True)

pca_button = st.button("Click button to run PCA analysis") # Give button a variable name
if pca_button: # Make button a condition.
    st.text("Start PCA analysis")
    st.write('PCA analysis with Leiden clustering')
    # KNN
    n_neighbors = 15 # Number of nearest neighbors for KNN graph
    knn_n_pcs = 50 # Number of principal components to use for finding nearest neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs)
    sc.pp.pca(adata, n_comps=50, use_highly_variable=False, svd_solver='arpack')
    sc.tl.leiden(adata)
    with rc_context({'figure.figsize': (10, 10)}):
        sc.pl.pca_scatter(adata,color=["leiden"],legend_loc='on data', legend_fontsize=10)
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
umap_button = st.button("Click button to run UMAP analysis") # Give button a variable name
if umap_button: # Make button a condition.
    st.text("Start UMAP analysis")
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
tsne_button = st.button("Click button to run t-SNE analysis") # Give button a variable name
if tsne_button: # Make button a condition.
    st.text("Start t-SNE analysis")
    sc.tl.tsne(adata, n_pcs=tsne_n_pcs)
    kmeans = KMeans(n_clusters=18, random_state=0).fit(adata.obsm['X_pca'])
    adata.obs['kmeans'] = kmeans.labels_.astype(str)
    with rc_context({'figure.figsize': (10, 10)}):
        st.text("t-SNE analysis result is here")
        fig = sc.pl.tsne(adata, color=["kmeans"],legend_loc='on data', legend_fontsize=10)
        st.pyplot()
st.session_state['adata'] = adata