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
import io
import matplotlib.pyplot as plt
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

image = Image.open('./images/deg-analysis.png')

st.image(image, use_column_width=True)

st.subheader('Differentially expressed gene analysis')


ranking_n_top_genes = 20 # Number of differential genes to compute for each cluster
adata = st.session_state['adata']
deg_button = st.button("Run DEG analysis") # Give button a variable name
if deg_button: # Make button a condition.
   st.text("Click button to run start DEG analysis")
   sc.tl.rank_genes_groups(adata, 'leiden', groups='all', reference='2', method='wilcoxon', key_added = "wilcoxon")
   st.text("DEG analysis result is here")
   sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, key="wilcoxon")
   img = io.BytesIO()
   plt.savefig(img, format='png')
   st.pyplot()
   st.download_button(label='Download DEG result (part 1)',
      data= img,
      file_name='deg_result_1.png',
      mime='image/png')
   sc.tl.rank_genes_groups(adata, groupby="leiden", n_genes=ranking_n_top_genes, groups='all', reference='rest', method='logreg')
   sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", groupby="leiden")
   img = io.BytesIO()
   plt.savefig(img, format='png')
   st.pyplot()
   st.download_button(label='Download DEG result (part 2)',
      data= img,
      file_name='deg_result_2.png',
      mime='image/png')
   sc.tl.rank_genes_groups(adata, 'majority_voting')
   sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key="wilcoxon", groupby="majority_voting", show_gene_labels=True)
   img = io.BytesIO()
   plt.savefig(img, format='png')
   st.pyplot()
   st.download_button(label='Download DEG result (part 3)',
      data= img,
      file_name='deg_result_3.png',
      mime='image/png')