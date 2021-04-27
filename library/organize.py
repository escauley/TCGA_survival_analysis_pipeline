# For manipulating file data.
import pandas as pd
import csv
import json
import sys
import requests

# For exiting the program.
import os

class organize:
    """ This class contains a set of functions for re-organizing gene expression and clinical data downloaded from the GDC data
    portal so that they can be used for downstream processing such as survival analysis.

    
    """

    def download_expression_data(self, sample_sheet, output_download):
        """ This function downloads expression data from GDC data portal specified in a downloaded sample sheet.

        First take a sample sheet downloaded from the GDC data portal with all samples to be downloaded.

        :param sample_sheet: The path to a sample sheet downloaded from the cart at the GDC data portal.
        :type sample_sheet: str
        :param output_directory: A path to the directory to populate with downloaded expression data.
        :type output_directory: str
        """

        # Set up a dictionary for collecting the informatiom form the sample sheet.
        collector = {}

        # Set up a counter for reading each row of the sample sheet file and skip the header.
        rowcount = 0

        # Open the sample sheet file.
        with open(sample_sheet) as sample_csv:
            csvreader = csv.reader(sample_csv, delimiter="\t", quotechar='"')
            for row in csvreader:
                rowcount += 1
                if rowcount == 1:
                    continue
                else:
                    # Save the file ID, project ID, and sample type.
                    fil_id = row[0]
                    proj_id = row[4]
                    sam_type = row[7].replace(" ", "-")

                    # Update collector dictionary
                    # Do the project id and sample type exist in the dictionary? If not, create it.
                    if not (proj_id in collector):
                        collector[proj_id] = {}
                    if not (sam_type in collector[proj_id]):
                        collector[proj_id][sam_type] = []

                    # Add the file Id to the collector dictionary.
                    collector[proj_id][sam_type].append(fil_id)

        # Download everything!
        data_endpt = "https://api.gdc.cancer.gov/data"

        # Set up the collector dictionary where the structure is project_id : sample_type : file_id.
        for proj_id in collector:

            # Create paths to establish destination folders.
            path1 = os.path.join(output_download, proj_id)
            if not (os.path.exists(path1)):
                os.mkdir(path1)
            for sam_type in collector[proj_id]:
                path2 = os.path.join(path1, sam_type)
                if not (os.path.exists(path2)):
                    os.mkdir(path2)

                # Make a list of all the files.
                fil_list = collector[proj_id][sam_type]

                # Create a file list of only unique names.
                uniq_fil_list = list(set(fil_list))

                params = {"ids": uniq_fil_list}

                response = requests.post(data_endpt, data=json.dumps(params),
                                         headers={"Content-Type": "application/json"})

                file_name = path2 + "/results.tar.gz"

                with open(file_name, "wb") as output_file:
                    output_file.write(response.content)

                os.mkdir(output_download + proj_id + '/' + 'logs/')

                log_file = output_download + proj_id + '/' + 'logs/' + proj_id + '.log'
                summary = str(proj_id + ',' + sam_type + ',' + str(len(fil_list)) + ',' + str(len(uniq_fil_list)))

                with open(log_file, 'w') as log:
                    log.write(summary)


    def map_ensg_to_genesymbol(self, ensg_mapping_file, uniprot_mapping_file, master_df):
        """ This function will map ENSG transcripts to gene symbol and Uniprot accessions.

        First the ENSG transcript is split into the base ENSG ID and transcript number,
        then creates a new column for gene symbol and maps the ENSG base ID to gene symbols using a provided mapping
        file. A column for Uniprot accession is then created and the gene symbol column is mapped using a provided
        mapping file so that every ENSG ID now has a corresponding gene symbol and Uniprot accession.

        :param ensg_mapping_file: The path to the file to map ENSG ID to Gene Symbol.
        :type ensg_mapping_file: str
        :param uniprot_mapping_file: The path to the file to map Uniprot accession to gene symbol
        :type uniprot_mapping_file: str
        :param master_df: A pandas dataframe with column for ensg_transcript
        :type master_df: A pandas dataframe object
        :return: A pandas dataframe with ENSG ID mapped
        :rtype: A pandas dataframe obj

        .. seealso:: For iterating over pandas dataframes:
            https://www.geeksforgeeks.org/iterating-over-rows-and-columns-in-pandas-dataframe/

        """

        # Load the ensg id mapping file.
        with open(ensg_mapping_file, "r") as ensg_file_handle:
            ensg_file_csv = csv.reader(ensg_file_handle, delimiter='\t')
            # Skip the header.
            next(ensg_file_csv)

            # Set up the mapping file dictionary.
            ensg_mapping_dict = {}

            # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
            for row in ensg_file_csv:
                ensg_mapping_dict[row[1]] = row[0]

        # Load the uniprot mapping file.
        with open(uniprot_mapping_file, "r") as uniprot_file_handle:
            uniprot_file_csv = csv.reader(uniprot_file_handle, quoting=csv.QUOTE_ALL)
            # Skip the header.
            next(uniprot_file_csv)

            # Set up the mapping file dictionary.
            uniprot_mapping_dict = {}

            # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
            for row in uniprot_file_csv:
                uniprot_mapping_dict[row[0]] = row[1]

        # Split the ENSG symbol column into a column for ENSG ID and ENSG Transcript
        master_df[['ensg_id', 'ensg_transcript_num']] = master_df.ensg_transcript.str.split('.', expand=True)

        print('Mapping to gene symbol')
        # Map the ENSG Id column to the gene symbols in the first mapping dict.
        master_df['gene_symbol'] = master_df['ensg_id'].map(ensg_mapping_dict)

        print('Mapping to uniprot accession')
        # Map the gene symbol column to the uniprot accession column in the second mapping dict.
        master_df['uniprotkb_ac'] = master_df['gene_symbol'].map(uniprot_mapping_dict)

        print('Mapping complete')

        return master_df

    def convert_fpkm_to_tpm(self, study_folder):
        '''

        :param study_folder: the folder containing samples to be processed
        :type study_folder: a path
        :return: 
        :rtype: a pandas dataframe
        '''


        # Arguments
        # ---------

        # study_folder: the folder containing samples to be processed

        # Returns
        # -------

        # A pandas dataframe with an additional column of TPM, converted from the FPKM column

        # Set up the path to the project directory.
        condition_directory = study_folder + 'Primary-Tumor/'

        # Go to directory of each project.
        os.chdir(condition_directory)

        # Set up a dictionary where keys are file ids and values are file names.
        file_id_dict = {}

        # Make a dictionary where keys are the file_id and values are the file names.
        with open('MANIFEST.txt', 'r') as manifest:

            # Skip the header.
            next(manifest)

            for row in manifest:
                split_row = row.split()

                # Save the file_id (first column)
                file_id = split_row[0]

                # Save the file name (second column)
                filename = split_row[1]
                if filename.endswith('.gz'):
                    filename = filename[:-3]

                # Story the file id and file name into the dictionary.
                file_id_dict[file_id] = filename

        # Define the columns to create in each sample dataframe.
        column_names = ['ensg_transcript', 'FPKM']

        # Go into each sample directory and calculate TPM for each gene.
        for folder, file in file_id_dict.items():
            print('processing sample: ' + folder)

            # Set up the path to the sample directory.
            sample_file = condition_directory + file

            # Load the sample file into a pandas dataframe.
            sample_df = pd.read_csv(file, sep='\t', names=column_names)

            # Sum the total FPKM for the sample.
            fpkm_sum = sample_df['FPKM'].sum()

            # Create a column for TPM and perform the calculation per row.
            sample_df['TPM'] = (sample_df['FPKM']/fpkm_sum)*1000000

            # Write the sample dataframe to the same sample directory.
            sample_df.to_csv(file, index=False)

            print(folder + ' complete')

    def write_out(self, final_dataframe, final_output_path):

        # Arguments
        # ---------

        # final_dataframe: the data frame to be written out.
        # final_output_path: the path to write the csv output file.

        # Returns
        # -------

        # A csv file from the pandas dataframe.

        print('Writing dataframe to ' + final_output_path)

        final_dataframe.to_csv(final_output_path, index=True)

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
                spl = line.split(sep=',')

                # Check that the log has four fields.
                if len(spl) != 4:
                    continue

                # Save the project name.
                proj_id = spl[0]

                # Add study and number of samples as key, value pairs to dictionary.
                if not (proj_id in study_list):
                    study_list.append(proj_id)

        # Go through each study and decompress the files.
        for proj_id in study_list:

            #Set up the path to the project directory.
            project_directory = data_folder + proj_id + '/Primary-Tumor/'

            # Go to directory of each project.
            os.chdir(project_directory)

            os.system('tar -zxvf results.tar.gz')

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

    def combine_tcga_readcounts(self, uncombined_input_folder):

        # Arguments
        # ---------

        # uncombined_input_folder: A full path to the directory that contains all the read count directories and filoes to combine.

        # Returns
        # -------

        # A pandas dataframe with all of the read count files combined and fields to denote file id and submitter id.

        # Designate the location of the files to be combined.
        os.chdir(uncombined_input_folder + 'Primary-Tumor/')

        # Set up a dictionary where keys are file ids and values are file names.
        file_id_dict = {}

        # Define the headers for the read count files and final master file.
        column_names = ['ensg_transcript', 'FPKM', 'TPM', 'file_id', 'file_name']

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
                if filename.endswith('.gz'):
                    filename = filename[:-3]

                # Story the file id and file name into the dictionary.
                file_id_dict[file_id] = filename

        print("finished making file id dictionary")

        # Use the file ID dictionary to edit and concatenate all read count files.
        for key, value in file_id_dict.items():

            # Check if the directory exists for that file ID.
            print('Checking for directory ' + key)
            if not os.path.isdir(key):
                # If it does not exist, move to the next file ID.
                continue

            # Load the read count file into a dataframe.
            df = pd.read_csv(value)

            # For all rows, make a constant value with the current file id and file name.
            df['file_id'] = key
            df['file_name'] = value

            print('Appending ' + value + ' to master df')

            # Add the read count dataframe to the master dataframe.
            master_df = master_df.append(df)

        return master_df

    def map_tcga_clinical_data(self, metadata, master_csv, clinical_folder):

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

        print('Mapping to metadata file')

        # Make a new column in the master df for the case id.
        # master_df.insert(3, 'case_id', [], True)

        # Map the metadata dictionary to the master dataframe.
        master_df['case_id'] = master_df['file_id'].map(metadata_dict)

        print('Metadata mapped')

        # Load the clinical data and gather the relavent data into mapping dictionaries.

        case_submitter_dict = {}
        age_dict = {}
        gender_dict = {}
        survival_days_dict = {}
        vital_status_dict = {}
        days_to_last_follow_up_dict = {}
        race_dict = {}
        ethnicity_dict = {}

        print('Mapping to clinical data')

        os.chdir(clinical_folder)

        # Uncompress clinical data file
        os.system('tar -xvzf *.tar.gz')

        # Load the clinical file.
        with open('clinical.tsv', "r") as clinical_handle:
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
                if row[0] not in days_to_last_follow_up_dict:
                    days_to_last_follow_up_dict[row[0]] = row[47]
                if row[0] not in race_dict:
                    race_dict[row[0]] = row[14]
                if row[0] not in ethnicity_dict:
                    ethnicity_dict[row[0]] = row[10]

        # Map the clinical dictionaries to the master dataframe.
        print('Mapping case_submitter_id')
        master_df['case_submitter_id'] = master_df['case_id'].map(case_submitter_dict)
        print('Mapping age_at_index')
        master_df['age_at_index'] = master_df['case_id'].map(age_dict)
        print('Mapping gender')
        master_df['gender'] = master_df['case_id'].map(gender_dict)
        print('Mapping vital_status')
        master_df['vital_status'] = master_df['case_id'].map(vital_status_dict)
        print('Mapping days_to_death')
        master_df['days_to_death'] = master_df['case_id'].map(survival_days_dict)
        print('Mapping days_to_last_follow_up')
        master_df['days_to_last_follow_up'] = master_df['case_id'].map(days_to_last_follow_up_dict)
        print('Mapping race')
        master_df['race'] = master_df['case_id'].map(race_dict)
        print('Mapping ethnicity')
        master_df['ethnicity'] = master_df['case_id'].map(ethnicity_dict)

        os.chdir('/mnt/c/Users/caule/PycharmProjects/survival_data')

        return master_df

    def generate_TPM_table(self, gene_list, cancer_list, data_folder):

        # Arguments
        # ---------

        # gene_list: a txt file of all of the genes to be included in the table
        # master_csv: a txt file of all of the TCGA cancer studies to be included in the table

        # Returns
        # -------

        # A pandas dataframe where rows are patient IDs, columns are genes, and values are TPM expression

        # Create lists for genes and cancers
        glist = []
        clist = []

        # Define the columns names for the master df.
        column_names = ['gene_symbol', 'TPM', 'tcga_patient_id', 'tcga_study']

        # Create the master pandas dataframe
        master_df = pd.DataFrame(columns=column_names)


        # Load the provided gene file into a list
        with open(gene_list, 'r') as gene_list_file:
            for row in gene_list_file:
                row = row.strip('\n')
                glist.append(row)


        # Load the provided cancer file into a list
        with open(cancer_list, 'r') as cancer_list_file:
            for row in cancer_list_file:
                row = row.strip('\n')
                clist.append(row)

        # Go through each cancer
        for cancer in clist:
            # Find the correct study folder and open the data file for that cancer
            os.chdir(data_folder + 'TCGA-' + cancer)

            # Define the cancer data file and open it as a Pandas dataframe.
            # Specified dtypes for each column to speed up reading large csv:
            # https://stackoverflow.com/questions/24251219/pandas-read-csv-low-memory-and-dtype-options
            cancer_data_file = 'human_cancer_TCGA_expression_survival_' + cancer + '.csv'
            cancer_data_df_full = pd.read_csv(cancer_data_file, quoting=csv.QUOTE_ALL,
                                              dtype={'uniprotkb_ac': 'string',
                                                     'gene_symbol': 'string',
                                                     'ensg_id': 'string',
                                                     'TPM': float,
                                                     'FPKM': float,
                                                     'tcga_patient_id': 'string',
                                                     'age': 'string',
                                                     'race': 'string',
                                                     'ethnicity': 'string'})
            # Create a column for the TCGA study.
            cancer_data_df_full['tcga_study'] = cancer

            # Take only the columns we will need.
            cancer_data_df_min = cancer_data_df_full[['gene_symbol', 'TPM', 'tcga_patient_id', 'tcga_study']]

            # Extract rows from the cancer data drame for the gene list
            cancer_data_df_filtered = cancer_data_df_min[cancer_data_df_min['gene_symbol'].isin(glist)]

            # Add the filtered rows to the master data frame.
            master_df = master_df.append(cancer_data_df_filtered, ignore_index=True)
            print('Summary of master df after processing ' + cancer)
            print(master_df.info())

        # Pivot the master dataframe so that each row is a patient ID and each column is a gene, values are TPM
        # Information about pivoting a table with Pandas:
        # https://stackoverflow.com/questions/47152691/how-to-pivot-a-dataframe
        print('Pivoting master dataframe')
        master_df = master_df.pivot_table(values='TPM', index=['tcga_patient_id', 'tcga_study'], columns='gene_symbol', aggfunc='mean')


        return master_df
