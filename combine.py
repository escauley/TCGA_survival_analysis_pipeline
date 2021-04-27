# To receive command line arguments.
import sys
import getopt

from library import organize
from library import FileUtils

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    #input_folder = ''
    #out_file_name = ''

    try:
        opts, args = getopt.getopt(argv, 'h:i:o:', ['input_folder=', 'out_file_name='])
    except getopt.GetoptError:
        print('combine_read_counts.py --input_folder [*ABSOLUTE PATH*] --out_file_name [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('combine_read_counts.py --input_folder [*ABSOLUTE PATH*] --out_file_name [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-i', '--input_folder'):
            input_folder = arg
        elif opt in ('-o', '--out_file_name'):
            out_file_name = arg

    # Instantiate organize.py
    org = organize.organize()

    # Create a new column for TPM.
    org.convert_fpkm_to_tpm(study_folder=input_folder)

    # Combine all read count files in a dataframe.
    combined_dataframe = org.combine_tcga_readcounts(uncombined_input_folder=input_folder)

    # Write combined dataframe to a master csv file.
    org.write_out(final_dataframe=combined_dataframe, final_output_path=out_file_name)

if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 combine.py -i /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-BRCA/ -o /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-BRCA/TCGA-BRCA_all_samples_TPM_1.csv