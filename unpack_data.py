# To receive command line arguments.
import sys
import getopt

from library import organize

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    #log_file = ''
    #data_folder = ''

    # To do: Change log file path so that it can take a TCGA study id instead.
    #   Combine with download data step

    try:
        opts, args = getopt.getopt(argv, 'h:l:d:', ['in_directory=', 'log_file_path=', 'data_folder_path='])
    except getopt.GetoptError:
        print('unpack_data.py --log_file_path [*ABSOLUTE PATH*] --data_folder_path [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('unpack_data.py --log_file_path [*ABSOLUTE PATH*] --data_folder_path [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-l', '--log_file_path'):
            log_file_path = arg
        elif opt in ('-d', '--data_folder_path'):
            data_folder_path = arg

    # Instantiate organize.py
    org = organize.organize()

    org.uncompress_tcga_hits(log_file=log_file_path, data_folder=data_folder_path)


if __name__ == '__main__':
    main(sys.argv[1:])

# Example command line call:
# python3 unpack_data.py -l /mnt/c/Users/caule/PycharmProjects/survival_data/logs/get_data_all_samples_BRCA.log -d /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/