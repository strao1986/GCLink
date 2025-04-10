# GCLink：Genetic Causality and Linking Molecular Mechanisms
GCLink is an integrated analytical pipeline designed to uncover the genetic causal relationships and shared molecular mechanisms between complex diseases, such as anxiety disorder (AD) and allergic rhinitis (AR). By employing large-scale epidemiological data, GWAS summary statistics, bulk-tissue and single-cell RNA sequencing eQTL datasets, GCLink provides a robust multi-stage approach for genetic causal inference and shared molecular mechanisms exploration.
## Overview
GCLink comprises four main phases (codes were provided with AD and AR as examples):
### 1.	Phenotypic Association Analysis
Employing large-scale epidemiological datasets (e.g., from the UK Biobank project) to establish logistic regression models for systematically evaluating phenotypic associations between complex diseases.
### 2.	Genetic Correlation and Causal Inference
(1) Employing Linkage Disequilibrium Score Regression (LDSC) program to assess genetic correlations between complex diseases.\
(2) Using MiXeR program to estimate polygenic overlap between complex diseases.\
(3) Performing bidirectional two-sample Mendelian Randomization (MR) analyses to investigate genetic causal correlations between complex diseases. By using the most suitable MR method for each MR analysis in different scenarios, we can minimize potential bias that caused by unsuitable methods. The most suitable MR method selection was described in detail in our previous studies (Rao _et al._, _J. Allergy Clin. Immunol._, 2025; Chen _et al._, _J. Transl. Med._, 2024; Rao _et al._, _Nutrients_, 2023).
### 3.	Shared Molecular Mechanisms Exploration
(1) Causal correlations of complex diseases risk with hematological traits and immune-related cells counts.  
(2) Exploring the shared genetic architecture between complex diseases at four levels:\
Genetic Loci: Identifying common risk SNPs and loci by the GWAS-PW program.\
Enriched Tissues: LDSC-SEG program for tissue-level shared heritability enrichments.\
Functional Genes: Summary-data-based Mendelian Randomization (SMR) framework for tissue-specific shared functional genes.\
Specific cell types: MAGMA.Celltyping program for tissue-specific shared risky cell types.
### 4.	Causal correlations of complex diseases risk with cognitive ability and behavioral symptoms changes
Assessing the genetic causal effects of complex diseases risk on changes of general cognitive ability and behavioral symptoms by bidirectional two-sample MR analyses, revealing the broader neuropsychological and behavioral impacts.
## Reference panels and preprocessed GWAS summary statistics downloading
The following resources are available for downloading from our Google Drive repository (Link: https://drive.google.com/drive/folders/1p3nH8z8tztblVUdQZzyPrZRj8l7Mtq9Y?usp=sharing):  
(1) Reference panels required in scripts running.\
(2) Preprocessed GWAS summary statistics for example diseases (AD and AR), hematological traits, immune-related cells counts, general cognitive ability and behavioral symptoms.
## Dependency
GCLink is primarily built using Python (≥3.7) and R (≥4.3.0). 
## Citation
If you utilize the GCLink pipeline-related scripts in your research, please cite the corresponding publications and related methodological references:  
(1)	Rao, S.T., J.W. Jiang, Y.Q. He, _et al._, Unraveling shared genetic architecture between multiple allergic diseases and anxiety disorder. _J. Allergy Clin. Immunol._, 2025. 155(2): p. AB229.\
(2)	Chen, X.T., S. Zhi, X.Y. Han, _et al._, A systematic two-sample and bidirectional MR process highlights a unidirectional genetic causal effect of allergic diseases on COVID-19 infection/severity. _J. Transl. Med._, 2024. 22(1): p. 94.

