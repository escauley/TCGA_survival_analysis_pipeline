# To receive command line arguments.
import sys
import getopt

from library import organize

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    # input_folder = ''
    # gene_list_input = ''
    # cancer_list_input = ''
    # output_file = ''

    try:
        opts, args = getopt.getopt(argv, 'h:i:g:c:o:', ['help', 'input_folder=', 'gene_list_input=', 'cancer_list_input=', '--ouput_file='])
    except getopt.GetoptError:
        print('generate_TPM_csv.py --help --input_folder [*ABSOLUTE PATH*] --gene_list_input [*ABSOLUTE PATH*] --cancer_list_input [*ABSOLUTE PATH*] --ouput_file [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('generate_TPM_csv.py --help --input_folder [*ABSOLUTE PATH*] --gene_list_input [*ABSOLUTE PATH*] --cancer_list_input [*ABSOLUTE PATH*] --ouput_file [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-i', '--input_folder'):
            input_folder = arg
        elif opt in ('-g', '--gene_list_input'):
            gene_list_input = arg
        elif opt in ('-c', '--cancer_list_input'):
            cancer_list_input = arg
        elif opt in ('-o', '--output_file'):
            output_file = arg

    # Instantiate organize.py
    org = organize.organize()

    # Create a new column for TPM.
    dataframe = org.generate_TPM_table(data_folder=input_folder, gene_list=gene_list_input, cancer_list=cancer_list_input)

    # Write combined dataframe to a master csv file.
    org.write_out(final_dataframe=dataframe, final_output_path=output_file)

if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 generate_TPM_csv.py -i /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/ -g /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/GT_gene_list.txt -c /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/GT_cancer_list.txt -o /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/GT_TPM_summary.csv