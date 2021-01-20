# To receive command line arguments.
import sys
import getopt

from library import organize

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    # meta_data_file = ''
    # master_csv_input = ''
    # clinical_tsv_input = ''
    # gene_symbol_mapping = ''
    # output_path = ''

    try:
        opts, args = getopt.getopt(argv, 'h:i:m:c:e:u:o:', ['master_csv_input=', 'metadata_file=', 'clinical_tsv_input=', 'ensg_mapping=', 'uniprot_mapping', 'output_path='])
    except getopt.GetoptError:
        print('map_to_metadata.py --metadata_file [*ABSOLUTE PATH*]' '--clinical_input_folder [*ABSOLUTE PATH*]' '--ensg_mapping [*ABSOLUTE PATH*]' '--uniprot_mapping [*ABSOLUTE PATH*]' '--output_path [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('map_to_metadata.py --metadata_file [*ABSOLUTE PATH*]' '--clinical_input_folder [*ABSOLUTE PATH*]' '--ensg_symbol_mapping [*ABSOLUTE PATH*]' '--uniprot_mapping [*ABSOLUTE PATH*]' '--output_path [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-i', '--master_csv_input'):
            master_csv_input = arg
        elif opt in ('-m', '--metadata_file'):
            metadata_file = arg
        elif opt in ('-c', '--clinical_input_folder'):
            clinical_input_folder = arg
        elif opt in ('-e', '--ensg_mapping'):
            ensg_mapping = arg
        elif opt in ('-u', '--uniprot_mapping'):
            uniprot_mapping = arg
        elif opt in ('-o', '--output_path'):
            output_path = arg

    # Instantiate organize.py
    org = organize.organize()

    master_df_with_clinical = org.map_tcga_clinical_data(master_csv=master_csv_input, metadata=metadata_file, clinical_folder=clinical_input_folder)

    master_df_with_gene_symbol = org.map_ensg_to_genesymbol(ensg_mapping_file=ensg_mapping, uniprot_mapping_file=uniprot_mapping, master_df=master_df_with_clinical)

    org.write_out(final_dataframe=master_df_with_gene_symbol, final_output_path=output_path)

if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 map_to_metadata.py
#   -i ../../OncoMX/survival_dataset/normalized_read_counts/TCGA-ESCA/TCGA-ESCA_all_samples_TPM.csv
#   -m ../../OncoMX/survival_dataset/normalized_read_counts/TCGA-ESCA/metadata.cart.2021-01-12.json
#   -c ../../OncoMX/survival_dataset/normalized_read_counts/TCGA-ESCA/clinical_info/
#   -e ../../OncoMX/survival_dataset/normalized_read_counts/ensgID_GeneSymbol_mapping.txt
#   -u ../../OncoMX/survival_dataset/normalized_read_counts/biomarker_masterlist.csv
#   -o ../../OncoMX/survival_dataset/normalized_read_counts/TCGA-ESCA/TCGA-ESCA_TPM_Survival.csv