# Single cell RNA-seq analysis using seurat

## Overview
This project performs an end-to-end single-cell RNA sequencing (scRNA-seq) analysis of human immune cells using the Seurat package in R.
The workflow includes quality control, normalization, dimensionality reduction, clustering, and marker gene identification.

## Dataset
The dataset was originally provided in `.h5ad` format and converted into 10X format using python for analysis in Seurat.

## Source:
https://cellxgene.cziscience.com/datasets

## Analysis Workflow
1. Load dataset in 10X format
2. Create Seurat object
3. Quality control filtering
4. Data normalization
5. Identification of highly variable genes
6. Data scaling
7. Principal Component Analysis (PCA)
8. Cell clustering
9. UMAP visualization
10. Marker gene identification

## Tools Used
- Python
- R
- Seurat
- dplyr

## Output Figures
- QC violin plots
- PCA elbow plot
- UMAP cluster visualization
- Marker gene heatmap
- Dot plot of marker genes

## Author
Mukesh Varathan
MSc Microbiology
