# For manipulating file data.
import pandas as pd
import csv
import json

# For exiting the program.
import sys
import os

# For the input file Regex.
import re

# For loading reading and loading files.
from . import FileUtils





class organize:

    # Class Description
    # -----------------

    # These methods are for taking the survival and expression datasets from TCGA and reorganizing their structure.

    def load_concat_files(self, data_file_list):

        # Arguments
        # ---------

        # data_file_list: a list of files to be joined together.

        # Returns
        # -------

        # A pandas data frame of the joined files.

        # Create a list to hold the pandas dataframes.
        df_all = []

        # Load each file and change column names
        for current_file in data_file_list:

            # Load the csv file into a pandas dataframe.
            df = pd.read_csv(current_file)

            # Rename the columns for the ENSG symbol, normal sample ID, and tumor sample ID.
            df.rename(columns={df.columns[0]: 'ensg_symbol', df.columns[1]: 'normal', df.columns[2]: 'tumor'}, inplace=True)

            # Take only the patient ID from the file name.
            patient_id = str(current_file)[-25:-13]

            print(patient_id)

            # Add a column for each row containing the patient ID and a column for mapped gene symbol.
            df['patient_id'] = patient_id

            # Split the ensg_symbol column into the ensg_id (for mapping) and the transcript number.
            df[['ensg_id', 'ensg_transcript_num']] = df.ensg_symbol.str.split('.', expand=True)

            # Add to the dataframe list.
            df_all.append(df)

        # Join all dataframes to one master dataframe.
        master_df = pd.concat(df_all, axis=0, ignore_index=True)

        return master_df

    def map_to_masterlist(self, mapping_file_list, data):

        # Arguments
        # ---------

        # data: a data frame of the data to be processed.

        # master_list_location: the path to the masterlist to use for mapping.

        # Returns
        # -------

        # A pandas dataframe of the data with a mapped column
        # A list of the IDs that were not mapped.

        # For iterating over pandas dataframes: https://www.geeksforgeeks.org/iterating-over-rows-and-columns-in-pandas-dataframe/

        for mapping_file in mapping_file_list:

            # Load the mapping file.
            with open(mapping_file, "r") as map_file_handle:
                map_file_csv = csv.reader(map_file_handle, delimiter='\t')
                # Skip the header.
                next(map_file_csv)

                # Set up the mapping file dictionary.
                mapping_dict = {}

                # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
                for row in map_file_csv:
                    mapping_dict[row[1]] = row[0]
                    #if row[1] in mapping_dict:
                        # Already there, additional term to add.
                        #mapping_dict[row[1]].append(row[0])
                    #else:
                        # Not present in dict yet so create it.
                        #mapping_dict[row[1]] = []
                        #mapping_dict[row[1]].append(row[0])

        # Map the unmapped data frame using the mapping dictionary.
        # Source: https://kanoki.org/2019/04/06/pandas-map-dictionary-values-with-dataframe-columns/
        data['gene_symbol'] = data['ensg_id'].map(mapping_dict)

        return data

    def transpose(self, untransposed_data):

        # Arguments
        # ---------

        # untransposed_data: a data frame of the data to be processed.

        # Returns
        # -------

        # A transposed dataframe with rows switched for columns.

        transposed_data = pd.transpose(data = untransposed_data)

        return transposed_data

    def write_out(self, final_dataframe, final_output_path):

        # Arguments
        # ---------

        # final_dataframe: the data frame to be written out.
        # final_output_path: the path to write the csv output file.

        # Returns
        # -------

        # A csv file from the pandas dataframe.

        final_dataframe.to_csv(final_output_path, index=False)

    def quote_csv(self, data_to_quote, output_name):

        # Arguments
        # ---------

        # data_to_quote: the data frame to be written out.

        # Returns
        # -------

        # A csv file with all values in double quotes.

        data_to_quote.to_csv('data/COAD/processed/COAD.csv', index=False)

    def uncompress_tcga_hits(self, log_file, data_folder):

        # Arguments
        # ---------

        # log_file: the list of studies downloaded as well as the number of samples per studied, generated from get_data_all_samples.py script.
        # data folder: the absolute path to the top level folder for the data to be uncompressed.

        # Returns
        # -------

        # None. Creates a folder with uncompressed sample files.

        # Create a list of the study names to process.
        study_list = []

        with open(log_file, 'r') as fil:
            # Go through each study and save the study name and number of samples.
            for line in fil:
                spl = line.split()

                # Check that the log has four fields.
                if len(spl) != 4:
                    continue

                # Save the project name.
                proj_id = spl[0]

                # Save the number of samples.
                num_sam = int(spl[2])

                # Add study and number of samples as key, value pairs to dictionary.
                if not (proj_id in study_list):
                    study_list.append(proj_id)

        # Go through each study and decompress the files.
        for proj_id in study_list:

            #Set up the path to the project directory.
            project_directory = data_folder + proj_id + '/Primary-Tumor/'

            # Go to directory of each project.
            os.chdir(project_directory)

            # Make a list of sample ids from the manifest.
            with open('MANIFEST.txt', 'r') as manifest:
                # Create a list to hold the sample ids.
                sample_list = []

                # Skip the header.
                next(manifest)

                # Got through the manifest and record each sample id.
                for line in manifest:

                    # The sample ids are in the first column in the manifest file.
                    manifest_spl = line.split()
                    sample_id = manifest_spl[0]

                    # Check if the sample id was already added to the list?
                    if sample_id not in sample_list:
                        sample_list.append(sample_id)
                    else:
                        continue

            # Go into each sample directory.
            for sample in sample_list:

                print('processing sample: ' + sample)

                # Set up the path to the sample directory.
                sample_directory = data_folder + proj_id + '/Primary-Tumor/' + sample

                # Go to the sample directory.
                os.chdir(sample_directory)
                
                # Uncompress the data.
                os.system('gunzip -k *.gz')

                print('sample complete')

    def combine_tcga_readcounts(self, uncombined_input_folder):

        # Arguments
        # ---------

        # uncombined_input_folder: A full path to the directory that contains all the read count directories and filoes to combine.

        # Returns
        # -------

        # A pandas dataframe with all of the read count files combined and fields to denote file id and submitter id.

        # Designate the location of the files to be combined.
        os.chdir(uncombined_input_folder)

        # Set up a dictionary where keys are file ids and values are file names.
        file_id_dict = {}

        # Define the headers for the read count files and final master file.
        column_names = ['ensg_transcript', 'FPKM', 'file_id', 'file_name']

        # Set up the master dataframe with header.
        master_df = pd.DataFrame(columns=column_names)

        # Load the MANIFEST into a list.
        with open('MANIFEST.txt', 'r') as manifest:

            next(manifest)

            for row in manifest:
                split_row = row.split()

                # Save the file_id (first column)
                file_id = split_row[0]

                # Save the file name (second column)
                filename = split_row[1]

                # Story the file id and file name into the dictionary.
                file_id_dict[file_id] = filename


        # Use the file id dictionary to edit and concatenate all read count files.
        for key, value in file_id_dict.items():

            # Load the read count file into a dataframe.
            df = pd.read_csv(value, sep='\t', names=column_names)

            # For all rows, make a constant value with the current file id and file name.
            df['file_id'] = key
            df['file_name'] = value

            # Add the read count dataframe to the master dataframe.
            master_df = master_df.append(df)

        return master_df

    def map_tcga_clinical_data(self, metadata, master_csv, clinical_tsv):

        # Arguments
        # ---------

        # metadata: the metadata file provided by TCGA.
        # master_csv: The csv of combined TCGA read counts.

        # Returns
        # -------

        # A pandas dataframe with clinical data from the metadata file mapped to every row.

        # Load the master data csv into a pandas dataframe.
        master_df = pd.read_csv(master_csv)

        # Create a dictionary to hold the metadata information.
        metadata_dict = {}

        # Load the metadata json file and relevant clinical data into the dictionary.
        with open(metadata) as f:
            metadata_json = json.load(f)

            # For each patient in the metadata file assign the file id to the case id in our dictionary.
            for patient in metadata_json:
                metadata_dict[patient['file_id']] = patient['associated_entities'][0]['case_id']

        # Make a new column in the master df for the case id.
        # master_df.insert(3, 'case_id', [], True)

        # Map the metadata dictionary to the master dataframe.
        master_df['case_id'] = master_df['file_id'].map(metadata_dict)

        # Load the clinical data and gather the relavent data into mapping dictionaries.

        case_submitter_dict = {}
        age_dict = {}
        gender_dict = {}
        survival_days_dict = {}
        vital_status_dict = {}
        race_dict = {}
        ethnicity_dict = {}


        # Load the clinical file.
        with open(clinical_tsv, "r") as clinical_handle:
            clinical_data = csv.reader(clinical_handle, delimiter='\t')
            # Skip the header.
            next(clinical_data)

            # Populate the mapping dictionaries with keys as case IDs and values as the clinical information.
            for row in clinical_data:
                if row[0] not in case_submitter_dict:
                    case_submitter_dict[row[0]] = row[1]
                if row[0] not in age_dict:
                    age_dict[row[0]] = row[3]
                if row[0] not in gender_dict:
                    gender_dict[row[0]] = row[11]
                if row[0] not in vital_status_dict:
                    vital_status_dict[row[0]] = row[15]
                if row[0] not in survival_days_dict:
                    survival_days_dict[row[0]] = row[9]
                if row[0] not in race_dict:
                    race_dict[row[0]] = row[14]
                if row[0] not in ethnicity_dict:
                    ethnicity_dict[row[0]] = row[10]

        # Map the clinical dictionaries to the master dataframe.
        master_df['case_submitter_id'] = master_df['case_id'].map(case_submitter_dict)
        master_df['age_at_index'] = master_df['case_id'].map(age_dict)
        master_df['gender'] = master_df['case_id'].map(gender_dict)
        master_df['vital_status'] = master_df['case_id'].map(vital_status_dict)
        master_df['vital_status'] = master_df['case_id'].map(vital_status_dict)
        master_df['days_to_death'] = master_df['case_id'].map(survival_days_dict)
        master_df['vital_status'] = master_df['case_id'].map(vital_status_dict)
        master_df['race'] = master_df['case_id'].map(race_dict)
        master_df['ethnicity'] = master_df['case_id'].map(ethnicity_dict)

        print(master_df.head())
        print(master_df.tail())






































