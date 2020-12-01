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

    try:
        opts, args = getopt.getopt(argv, 'h:i:m:c:', ['master_csv_input=', 'metadata_file=', 'clinical_tsv_input='])
    except getopt.GetoptError:
        print('map_to_metadata.py --metadata_file [*ABSOLUTE PATH*]' '--clinical_tsv_input [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('map_to_metadata.py --metadata_file [*ABSOLUTE PATH*]' '--clinical_tsv_input [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-i', '--master_csv_input'):
            master_csv_input = arg
        elif opt in ('-m', '--metadata_file'):
            metadata_file = arg
        elif opt in ('-c', '--clinical_tsv_input'):
            clinical_tsv_input = arg

    # Instantiate organize.py
    org = organize.organize()

    org.map_tcga_clinical_data(master_csv=master_csv_input, metadata=metadata_file, clinical_tsv=clinical_tsv_input)


if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 map_to_metadata.py -i /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-PRAD/TCGA-PRAD_all_samples_FPKM.csv -m /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-PRAD/metadata.cart.2020-11-25.json -c /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-PRAD/clinical_info/clinical.tsv