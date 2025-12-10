# PTB_BMI_protein_study
Code and analysis pipeline for the manuscript:
**Protein Networks are Influenced by Maternal BMI and Differentiate Preterm Birth Types (Lopez Zapana et al., 2025)**

## Overview
This repository contains the analysis scripts, helper functions, and figure-generation code used in the manuscript. The study integrates SomaScan v4.1 aptamer-based proteomics with:

- Principal component analysis (PCA)
- Elastic Net–regularized partial least squares discriminant analysis (PLSDA)
- Reactome pathway enrichment
- Graphical Lasso–based partial correlation network modeling
- In-silico perturbation analysis

The goals of the analysis are to:
1. Identify protein signatures distinguishing spontaneous preterm birth (sPTB), medically-indicated preterm birth (mPTB), and term deliveries.
2. Characterize how maternal early-pregnancy BMI modifies the proteomic architecture of each PTB subtype.
3. Identify proteins that act as key intermediates linking BMI to PTB risk.

All analyses were performed in R (version 4.3.2 or later).

## Repository Structure
```
PTB_BMI_protein_study/
│
├── Figure_1/               PCA and metadata associations
├── Figure_2/               Three-class PLSDA and multivariate loadings
├── Figure_3/               Subtype-specific PLSDA and Reactome enrichment
├── Figure_4/               BMI differential expression and BMI×PTB interaction
├── Figure_5/               Graphical Lasso partial correlation networks
├── Figure_6/               Perturbation (node removal) analysis
│
├── Sup_Figs/               Supplementary figure scripts
│
├── helpful_functions/      Helper R functions sourced by all scripts
│       data_loading.R
│       plotting_utils.R
│       plsda_modeling.R
│       pathway_analysis.R
│       network_glasso.R
│       perturbation_analysis.R
│
├── .gitignore
└── README.md
```

## Data Availability
Raw SomaScan proteomic data will be deposited in ImmPort under accession SDY3306.  
Intermediate processed files used to generate manuscript figures are included where allowed or automatically produced by scripts.
Other data is available upon request. 

## Analytical Workflow Summary
### 1. PCA
PCA on 6,926 log2-transformed proteins; associations with maternal metadata (BMI, PTB status, GA at draw).

### 2. Elastic Net + PLSDA
Elastic Net repeated 100×; three‑class and subtype‑specific PLSDA models; empirical null validation.

### 3. Reactome Pathway Enrichment
Proteins correlated with latent variables tested for pathway enrichment.

### 4. BMI Differential Expression & Interaction (limma)
BMI <30 vs ≥30 comparisons and BMI×PTB interaction assessed using likelihood‑ratio tests.

### 5. Graphical Lasso Partial Correlation Networks
Stability‑selected sparse networks identifying direct BMI→protein→PTB associations.

### 6. Perturbation Analysis
Node removal from precision matrix to identify key intermediates influencing BMI–PTB correlation.

## Outputs
Running all figure scripts produces:
- High‑resolution PDFs
- Differential expression tables
- Partial correlation matrices
- Perturbation results

## Citation
If you use this code, please cite:
Lopez Zapana PA, DeBolt CA, et al. Protein Networks are Influenced by Maternal BMI and Differentiate Preterm Birth Types. 2025.

## Contact
For questions, please open a GitHub Issue or contact the corresponding author.
