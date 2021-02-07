# To receive command line arguments.
import sys
import getopt

from library import organize

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    #sample_sheet_file = ''
    #download_folder = ''

    try:
        opts, args = getopt.getopt(argv, 'h:s:d:', ['sample_sheet_file=', 'download_folder='])
    except getopt.GetoptError:
        print('unpack_data.py --sample_sheet_file [*ABSOLUTE PATH*] --download_folder [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('unpack_data.py --sample_sheet_file [*ABSOLUTE PATH*] --download_folder [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-s', '--log_file_path'):
            sample_sheet_file = arg
        elif opt in ('-d', '--download_folder'):
            download_folder = arg

    # Instantiate organize.py
    org = organize.organize()

    org.download_expression_data(sample_sheet=sample_sheet_file, output_download=download_folder)


if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 download_data.py -s /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-UVM/gdc_sample_sheet.2021-02-07.tsv -d /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-UVM/