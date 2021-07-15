# import datetime
import sys
import logging
import os
logger = logging.getLogger('main')

mapitope = {
    'R': 'R',
    'K': 'R',
    'E': 'E',
    'D': 'E',
    'S': 'S',
    'T': 'S',
    'I': 'I',
    'L': 'I',
    'V': 'I',
    'Q': 'Q',
    'N': 'Q',
    'Y': 'Y',
    'F': 'Y',
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'H': 'H',
    'M': 'M',
    'P': 'P',
    'W': 'W'
}

def convert_to_mapitope(sequence):
    return ''.join([mapitope[c] for c in sequence])


def convert_faa_file_by_maptope(input_file_path, output_file_path, done_file_path, argv='no_argv'):
    with open(input_file_path,'r') as input_file, open(output_file_path,'w+') as output_file:
        input_lines = input_file.readlines()
        for line in input_lines:
            if line[0] == '>':
                output_file.write(line)  
            else:
                output_file.write(convert_to_mapitope(line[:-1])+'\n')

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='path to reads_filtration folder')
    parser.add_argument('output_path', help='path to barcode2sample file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    convert_faa_file_by_maptope(args.input_path, args.output_path, args.done_file_path,sys.argv)
    