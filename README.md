# scRNAseq analysis GUI using R

**GUI Overview:**
This project provides an interactive R Shiny GUI for analyzing single-cell RNA-seq data. The interface is organized into modular components, making it easy to follow the standard scRNA-seq analysis workflow:

* Upload Module – Supports input of raw count data in 10X Genomics (folder format) or H5 files. Users can add multiple samples and assign custom names.

* QC & Filtering Module – Allows quality control by filtering low-quality cells/genes using customizable thresholds.

* Normalization & Scaling – Implements different normalization strategies (e.g., Log, CLR, RC) and optional scaling for downstream analysis.

* Dimensionality Reduction – Performs PCA and visualizes results, with tools like the Elbow Plot to guide PC selection.

* Clustering - Perform clustering of cells to identify distinct populations and visualize them in reduced dimensions (e.g., PCA/UMAP).

* Cell Annotation - Annotate clusters based on known biological information or reference datasets to assign meaningful cell types.

* Gene Highlighting - Highlight and visualize the expression of specific genes across clusters and conditions.

* Differentially Expressed Genes (DEGs) - Identify and explore marker genes by running DEG analysis between selected groups of cells.

## How to run the program

Before running the program, ensure that you have:
* R environment installed.
* The following R libraries installed:
  ```shiny, shinyFiles, Seurat, SingleCellExperiment, SingleR, celldex, dplyr, ggplot2, biomaRt, shinycssloaders, hdf5r, colourpicker.```
<br><br>

For your convevience, run the following command in your R evironment to download the required libraries:<br>
```
install.packages(setdiff(
  c("shiny", "shinyFiles", "Seurat", "SingleCellExperiment", "SingleR", "celldex", "dplyr", "ggplot2", "biomaRt", "shinycssloaders", "hdf5r", "colourpicker"),
  rownames(installed.packages())
))
```


Fisrt, clone or download the repository to your PC.<br>
Next, 
**running instructions

Once you launch the app, the following screen should appear in your browser:

<img width="600" height="300" alt="image" src="https://github.com/user-attachments/assets/50e59d2c-fe71-4fde-948c-a143f345eda3" />


<br>Now, you can proceed with your analysis according to the data type.

## Uploading data

**The GUI supports scRNA-seq data in three formats: .rds, 10X, and .h5. <br>**
After selecting the data type, provide the path to your dataset. <br><br>
The path must follow this format: ```your\path\to\data.rds```<br><br>
**Important:** Make sure the file path contains no spaces or parentheses.<br>
You can easily copy the correct path by right-clicking the file and selecting "Copy as path" (Windows).


## How to process .rds files

In the main upload screen, select .rds data type, insert path to data and press ```Load Data```.<br>
The main analysis page would appear:<br><br>
<img width="800" height="500" alt="image" src="https://github.com/user-attachments/assets/53f016c5-b8d9-436f-b646-1c32c546b32f" />
<br>

Use the navbar at the top of the page to navigate between the different analysis pages:
* Cell annotation - Possible for Human or mouse samples
* Gene highlighting - Main analysis page
* Differentially expressed genes
* rds Download - Download updated version of .rds file


## How to process 10X or H5 files
