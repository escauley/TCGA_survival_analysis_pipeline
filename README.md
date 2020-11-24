# TCGA_survival_analysis_pipeline

## Summary
A series of shell and python scripts with options that can be called from the command to download read count data from the GDC portal and organizing the data so that it is mapped to patient clinical information for survival analysis.


## Pipeline steps

### Download

Top level script: get_data_all_samples.sh

Library Script: get_data_all_samples.py

Summary: Adapted from the Bioxpress pipeline Downloader step. Hard coded to download normalized read count data from TCGA studies depending on the sample sheet provided. 


### Uncompress files

Top level script: unpack_data.py

Options
-h or --help 
  - Display the options on the command line. 
  
-l or --log_file_path
  - The full path to the log file generated in the Download step
  
-d or --data_folder_path
  - The full path to the folder containing the read count files downloaded in the Download step
  
Library Script: organize.py
- Uses the function uncompress_tcga_hits




