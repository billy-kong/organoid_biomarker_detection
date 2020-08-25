# organoid_biomarker_detection
# Network-based biomarker detection from organoid models to predict drug response in cancer patients
Source code to reproduce the paper "Network-based machine-learning in bladder and colorectal organoid models predicts anti-cancer drug efficacy in patients", Kong et al

# Requirements
- pandas
- matplotlib
- numpy
- scipy
- sklearn
- lifelines
- gseapy


# Code (for python)
- "run_ssGSEA.py" to generate pathway level expression profiles using single sample GSEA (ssGSEA) tool (gseapy)
- "single_pathway_prediction.py" to predict drug response in cancer patients using a single pathway
- "multiple_pathway_prediction.py" to predict drug response in cancer patients using multiple pathways

# Demo
Code for drug response prediction of 5fluorouracil-treated colorectal cancer patients using colorectal cancer organoids


# Network proximity was calculated using codes from
- 'Uncovering disease-disease relationships through the incomplete interactome' Menche et al, Science, 2015
- 'Network-based in silico drug efficacy screening' Emre et al, Nature Communications, 2016
- https://github.com/emreg00/toolbox

