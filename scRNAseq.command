#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
cd "scRNAseq_GUI-main/scRNA_app"
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
