# Checker-Toxicity

### Overview
This repository contains the custom toxicity checks for bight

All of the code stored here is specific to bight toxicity checks only. The code is intended to be run standalone. The only packages you would need are what is listed in the import section at the top of each file (plus xlsxwriter). When the code is finished it is imported into the main checker application (SCCWRP/Checker)

### Instructions and Information
Put the data file in the data folder, and run from the command line:  
> python3 main.py -f /path/to/toxdata.xlsx

The report will go to output/report.xlsx (if the output folder does not exist, the script will create it). This report will contain the original dataframes, result, batch and water quality (wq) along with their respective error reports.
This report.xlsx file will also have the toxicity summary report generated from the raw data

There is an example data file in the data folder, called example.xlsx, the correct output for that data file is in validation/example_report.xlsx

Be sure that your input data has the same column headers as example.xlsx, and that the data is formatted in the same fashion

To print brief instructions, you can type 
> python3 main.py -h


