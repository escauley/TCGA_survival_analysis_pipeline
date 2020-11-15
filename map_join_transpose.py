# To receive command line arguments.
import sys
import getopt

from library import organize

def main(argv):

    # Get the command line arguments.
    # The parameters we're looking for.
    #mapping_file = ''
    #in_regex = ''
    #in_directory = ''

    try:
        opts, args = getopt.getopt(argv, 'h:i:r:m', ['in_directory=', 'in_regex=', 'mapping_file='])
    except getopt.GetoptError:
        print ('map_join_transpose.py --in_directory [*ABSOLUTE PATH*] --in_regex [Any Regex] --mapping_file [*ABSOLUTE PATH*]')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print ('map_join_transpose.py --in_directory [*ABSOLUTE PATH*] --in_regex [Any Regex] --mapping_file [*ABSOLUTE PATH*]')
            sys.exit()
        elif opt in ('-i', '--in_directory'):
            in_directory = arg
        elif opt in ('-r', '--in_regex'):
            in_regex = arg
        elif opt in ('-m', '--mapping_file'):
            mapping_file = arg

    # Instantiate organize.py
    org = organize.organize()

    # Load the files provided using the regex.

    files = org.read_files(input_directory=in_directory, regex=in_regex)

    bcos = cts.map_to_masterlist(file_locations=files)

    # Load the schema.

    schema = fileutils.read_files(input_directory=schema_location, regex=schema_name)
    print(schema)

    with open(schema[0]) as f:
        schema = json.load(f)

    cts.create_comparison_file(p_bcos=bcos, incoming_schema=schema)

    #for key, value in bcos.items():
        #print(key)
        #print(value)
        #print('++++++++++++++++++++++++++++++++++++\n\n\n\n\n\n')


if __name__ == '__main__':
    main(sys.argv[1:])

# python3 compare.py -i schema_convert_test_files/ -r 'bco_set_*.txt' -s schemas/ -n '2791object.json'