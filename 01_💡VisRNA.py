######################
# Import libraries
######################

import streamlit as st
from PIL import Image

######################
# Page Title
######################
st.write("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Arimo');
html, body, [class*="css"]  {
   font-family: 'Arimo';
}
</style>
""", unsafe_allow_html=True)
image = Image.open('./images/main-page.png')

st.image(image, use_column_width=True)

image = Image.open('./images/Workflow.png')
# st.write("## Visualization of tissue heterogeneity of non-small cell lung cancer by single-cell transcriptomics (VisRNA)")
st.image(image, use_column_width=True)
st.markdown('<div style="text-align: justify;">Lung cancer is one of the most significant causes\
   of death among all diseases. Specifically, non-small cell lung cancer (NSCLC) accounts for 87% of this daunting disease.\
   Single-cell RNA-sequencing (scRNA-seq) measures the transcriptomics profiles for each cell,\
   which facilitate the understanding of disrupted genes and mechanisms at the single-cell level.\
   This research project developed an efficient visualization platform (VisRNA) to conduct statistical\
   and functional analysis using scRNA-seq data to investigate the essential driver genes for understanding\
   the mechanisms of NSCLC. VisRNA enables interactive data uploading, scRNA-seq data processing and visualization,\
   and functional analysis. The VisRNA platform first performed multiple dimensionality reduction techniques of\
   scRNA-seq data, followed by automatic cell type annotation by machine learning algorithms. Differentially\
   expressed genes (DEGs) analysis and gene enrichment analysis were performed. Fifteen cell types were annotated\
   with DEGs when applying VisRNA to analyze NSCLC scRNA-seq data. In addition, several key driver genes, such as\
   RPL32 and DRAM1, were identified to be associated with NSCLC. The VisRNA platform provides an efficient and\
   user-friendly platform to analyze scRNA-seq data, which will promote cancer biomarker discovery scRNA-seq-based\
   biomedical research.</div>', unsafe_allow_html=True)
st.markdown('<p>&nbsp</p>', unsafe_allow_html=True)
st.markdown('<div style="text-align: justify;">Data was obtained from the CancerSEM database.</div>', unsafe_allow_html=True)