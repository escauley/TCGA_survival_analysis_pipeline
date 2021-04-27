# TCGA_survival_analysis_pipeline

## Summary
A series of shell and python scripts with options that can be called from the command to download read count data from the GDC portal and organizing the data so that it is mapped to patient clinical information for survival analysis.


## Pipeline steps


### Download

#### Top level script
download_data.py

Command Line Options
-h or --help
  - Display the options on the command line

-d or --download_folder
  - The main folder to download the read count data (a folder for the specific TCGA study will be generated in this folder)

#### Library Script 
organize.py

#### Summary
Adapted from the Bioxpress pipeline Downloader step. 

#### To Run
Go to the GDC data portal (https://portal.gdc.cancer.gov/) and add all the samples you wish to download to the cart. 

To add samples to the cart:
- Click the Advanced Search Tab and enter:

cases.project.program.name in ["TCGA"] and cases.project.project_id in ["TCGA-BLCA"] and files.analysis.workflow_type in ["HTSeq - FPKM"] and files.data_type in ["Gene Expression Quantification"] and files.experimental_strategy in ["RNA-Seq"] and cases.samples.sample_type in ["Primary Tumor"]

Download the sample sheet, the manifest, the metadata file, and the clinical data (as a tsv file) from the cart and onto your machine where you you will know the full path for referencing these files in later scripts. 

Run the download_data.py script
Example usage:
python download_data.py -d /data/

***Note***
After the download completes, move the sample sheet, manifest, metadata file, and clinical data into the TCGA study folder (example: TCGA-BLCA).

Make a folder called "clinical_info" and put the clinical data file withint the clinical info folder.

### Uncompress files

#### Top level script: 
unpack_data.py

Command Line Options
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


***Note***
If the download step did not generate a file named "MANIFEST.txt" in the same folder as the zipped data, this step will encounter an error and not run. In this case, make a copy of the manifest file downloaded in the previous step and name this copy MANIFEST.txt, then move this file to the same folder as the zipped read count data.
 - An example path for this file is /data/TCGA-BLCA/Primary-tumor/MANIFEST.tsxt

### Combine files

#### Note: 
In this step, FPKM values are also converted to TPM values for each sample using the formula:

TPM = (FPKM gene / Sum of FPKMs for all genes in the sample) * 1,000,000

TPM values are added to the original, uncompressed sample files before combining all samples into the master file. 

#### Top level script: 
combine.py

Command Line Options
-h or --help 
  - Display the options on the command line. 
  
-i or --input folder
  - The full path to the top level folder containing all data in the study being processed
  
-o or --out_file_name
  - The full path and name of the out file containing all read count/expression data
  
#### Library Script: 
organize.py
- Uses the function convert fpkm_to_tpm, combine_tcga_readcounts, and write_out

#### Example command line usage: 

python3 combine.py -i /TCGA-PRAD/Primary-Tumor/ -o TCGA_PRAD/all_samples_combined.csv


### Map to metadata, gene symbols, and clinical data

#### Top level script: 
map_to_metadata.py

Note: This step is used to map the sample id contained in expression files to the patient id found in the metadata file, then use the patient id from the metadata file to map clinical data in the clinical data file

Command Line Options
-h or --help 
  - Display the options on the command line. 
  
-i or --master_csv_input
  - The full path to the file containing all samples merged together (the table output from the previous Combine step)
  
-m or --metadata_file
  - The full path to the metadata file downloaded from the GDC data portal
  
-c or --clinical_tsv_output
  - the full path to the clinical tsv file downloaded from the GDC data portal
  
-g or --gene_symbol_mapping
  - The full path to the mapping file with ensg IDs and corresponding gene symbols, found in the mapping folder in the github repository
  
-o or --output_path
  - The full path and file name for the output file with all rows mapped 
  
#### Library Script: 
organize.py
- Uses the functions map_tcga_clinical_data, map_ensg_to_genesymbol, and write_out

#### Example command line usage: 

python3 map_to_metadata.py -i TCGA_PRAD/all_samples_combined.csv -m TCGA-PRAD/metadata.cart.2020-11-25.json -c TCGA-PRAD/clinical.tsv -g mapping/ensgID_GeneSymbol_mapping.txt -o TCGA-PRAD/TCGA-PRAD_FPKM_Survival.csv

