# For manipulating file data.
import pandas as pd
import csv

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

        # unmapped_data: a data frame of the data to be processed.

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

    def write_out(self, final_dataframe):

        # Arguments
        # ---------

        # final_dataframe: the data frame to be written out.
        # output_path: the path to write the csv output file.

        # Returns
        # -------

        # A csv file from the pandas dataframe.

        final_dataframe.to_csv('data/COAD/processed/COAD.csv', index=False)

    def quote_csv(self, data_to_quote, output_name):

        # Arguments
        # ---------

        # data_to_quote: the data frame to be written out.

        # Returns
        # -------

        # A csv file with all values in double quotes.

        data_to_quote.to_csv('data/COAD/processed/COAD.csv', index=False)





