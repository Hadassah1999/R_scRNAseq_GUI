@echo off
REM Get the folder where this script is located
SET SCRIPT_DIR=%~dp0

REM Change to that folder
cd /d "%SCRIPT_DIR%"

REM Change into scRNA_app subfolder
cd scRNA_app

REM Run the Shiny app with R
Rscript -e "shiny::runApp('.', launch.browser = TRUE)"

pause
