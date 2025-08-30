# scRNAseq analysis GUI using R
GUI Overview

* This project provides an interactive R Shiny GUI for analyzing single-cell RNA-seq data. The interface is organized into modular components, making it easy to follow the standard scRNA-seq analysis workflow:

* Upload Module – Supports input of raw count data in 10X Genomics (folder format) or H5 files. Users can add multiple samples and assign custom names.

* QC & Filtering Module – Allows quality control by filtering low-quality cells/genes using customizable thresholds.

* Normalization & Scaling – Implements different normalization strategies (e.g., Log, CLR, RC) and optional scaling for downstream analysis.

* Dimensionality Reduction – Performs PCA and visualizes results, with tools like the Elbow Plot to guide PC selection.

* Clustering - Perform clustering of cells to identify distinct populations and visualize them in reduced dimensions (e.g., PCA/UMAP).

* Cell Annotation - Annotate clusters based on known biological information or reference datasets to assign meaningful cell types.

* Gene Highlighting - Highlight and visualize the expression of specific genes across clusters and conditions.

* Differentially Expressed Genes (DEGs) - Identify and explore marker genes by running DEG analysis between selected groups of cells.

## How to run the program
1) first, clone repository


## How to process .rds files

## How to process 10X or H5 files
