######################
# Import libraries
######################

import streamlit as st
from PIL import Image

######################
# Page Title
######################

image = Image.open('./images/Workflow.png')

st.image(image, use_column_width=True)

st.write("""
## Visualization of tissue heterogeneity of non-small cell lung cancer by single-cell transcriptomics (VisRNA)

Lung cancer is one of the most significant causes of death among all diseases. 
Specifically, non-small cell lung cancer (NSCLC) accounts for 87% of this daunting disease. 
VisRNA enables interactive data uploading, scRNA-seq data processing and visualization, 
and functional analysis. Single-cell RNA-sequencing (scRNA-seq) is a technique that measures 
the transcriptomics profiles for each cell. This research aims to develop an efficient visualization platform (VisRNA) 
to facilitate downstream analysis using scRNA-seq data, which is critical to find essential driver genes for 
understanding the mechanisms of NSCLC. The VisRNA platform includes three core components, dimension reduction of high-dimensional 
scRNA-seq data, cell type annotation by known marker genes and machine learning, and functional annotation of differentially expressed 
genes (DEGs) for NSCLC. Fifteen cell types were annotated with DEGs when applying VisRNA to analyze NSCLC scRNA-seq. 
In addition, several key driver genes were identified to be associated with NSCLC. The visRNA platform will promote biomarker discovery based on scRNA-seq in biomedical research.


Data was obtained from the CancerSEM database.
""")