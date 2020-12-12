# TCGA_survival_analysis_pipeline

## Summary
A series of shell and python scripts with options that can be called from the command to download read count data from the GDC portal and organizing the data so that it is mapped to patient clinical information for survival analysis.


## Pipeline steps

### Download

Top level script: get_data_all_samples.sh

Library Script: get_data_all_samples.py

Summary: Adapted from the Bioxpress pipeline Downloader step. Hard coded to download normalized read count data from TCGA studies depending on the sample sheet provided. This will be updated to take command line arguments in the future updates.

### To Run
Go to the GDC data portal (https://portal.gdc.cancer.gov/) and add all the samples you wish to download to the cart. 

#### Example: 
Click the Advanced Search Tab and enter

files.analysis.workflow_type in ["FPKM"]  and files.data_type in ["Gene Expression Quantification"] and cases.samples.sample_type in ["Primary Tumor"] and cases.project.program.name in ["TCGA"] and cases.project.study.name in ["PRAD"]

Download the sample sheet, the manifest, and the clinical data (as a csv file) from the cart. 

Open get_data_all_samples.sh andinput the full path to the sample sheet downloaded from the GDC data portal (the second argument called to the python script).

Run the shell script get_data_all_samples.sh

### Uncompress files

#### Top level script: 
unpack_data.py

#### Options
-h or --help 
  - Display the options on the command line. 
  
-l or --log_file_path
  - The full path to the log file generated in the Download step
  
-d or --data_folder_path
  - The full path to the folder containing the read count files downloaded in the Download step
  
#### Library Script: 
organize.py
- Uses the function uncompress_tcga_hits

#### Example command line usage: 

python3 unpack_data.py -l logs/get_data_all_samples.log -d TCGA_data/normalized_read_counts/




