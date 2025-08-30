#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"
cd "scRNA_app"
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"
