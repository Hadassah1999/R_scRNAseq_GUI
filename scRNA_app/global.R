# List of required libraries
libs <- c(
  "shiny",
  "shinyFiles",
  "Seurat",
  "SingleCellExperiment",
  "SingleR",
  "celldex",
  "dplyr",
  "ggplot2",
  "biomaRt",
  "shinycssloaders",
  "hdf5r",
  "colourpicker"
)

# Install any that are missing
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib)
  }
}

# Load the libraries
library(shiny)
library(shinyFiles)
library(Seurat)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(shinycssloaders)
library(hdf5r)
library(colourpicker)
