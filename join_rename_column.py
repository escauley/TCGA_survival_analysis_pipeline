# To receive command line arguments.
import sys
import getopt

from library import organize
from library import FileUtils

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for:
    #in_regex = ''
    #in_directory = ''
    #mapping_file_directory = ''
    #mapping_file_regex = ''

    try:
        opts, args = getopt.getopt(argv, 'h:i:r:d:m:', ['in_directory=', 'in_regex=', 'mapping_file_directory=', 'mapping_file_regex='])
    except getopt.GetoptError:
        print ('join_rename_column.py --in_directory [*ABSOLUTE PATH*] --in_regex [Any Regex] --mapping_file_directory [*ABSOLUTE PATH*] --mapping_file_regex [Any Regex]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print ('join_rename_column.py --in_directory [*ABSOLUTE PATH*] --in_regex [Any Regex] --mapping_file_directory [*ABSOLUTE PATH*] --mapping_file_regex [Any Regex]')
            sys.exit()
        elif opt in ('-i', '--in_directory'):
            in_directory = arg
        elif opt in ('-r', '--in_regex'):
            in_regex = arg
        elif opt in ('-d', '--mapping_file_directory'):
            mapping_file_directory = arg
        elif opt in ('-m', '--mapping_file_regex'):
            mapping_file_regex = arg

    # Instantiate organize.py
    org = organize.organize()

    # Instantiate FileUtils
    fu = FileUtils.FileUtils()

    # Load the files provided using the regex.
    map_files = fu.read_files(input_directory=mapping_file_directory, regex=mapping_file_regex)
    files = fu.read_files(input_directory=in_directory, regex=in_regex)

    # Join the files.
    unmapped_master_file = org.load_concat_files(data_file_list=files)

    # Map the master_file.
    mapped_master_file = org.map_to_masterlist(mapping_file_list=map_files, data=unmapped_master_file)

    # Write out the master file.
    org.write_out(final_dataframe=mapped_master_file)


if __name__ == '__main__':
    main(sys.argv[1:])

# python3 join_rename_column.py -i data/COAD/ -r '*.csv' -d mapping/ -m '*.txt'