import subprocess
import logging
import datetime
import os
import sys
from time import time
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, measure_time


def cluster_sequences(fasta_file, output_prefix, done_file_path, threshold, word_length,
                      throw_sequences_shorter_than, cluster_algorithm_mode, argv='no_argv'):

    verify_file_is_not_empty(fasta_file)

    logger.info(f'{datetime.datetime.now()}: clustering sequences in {fasta_file}')

    # TODO: module load CD hit
    cmd = f'cd-hit -i {fasta_file} ' \
          f'-o {output_prefix} ' \
          f'-c {threshold} ' \
          f'-n {word_length} ' \
          f'-l {throw_sequences_shorter_than} ' \
          f'-g {cluster_algorithm_mode}'
    logger.info(f'Starting CD-hit. Executed command is:\n{cmd}')
    subprocess.call(cmd, shell=True)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(f'{output_prefix}.clstr')

    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')



if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', help='A fasta file to collapse for unique sequences and their counts')
    parser.add_argument('output_prefix', help='A file prefix in which the 2 result files will use as an output path')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('sample_name', type=str, help='The name of the sample that cluster on')
    
    parser.add_argument('--threshold', default='0.5', help='Minimal sequence similarity threshold required',
                        type=lambda x: float(x) if 0.4 <= float(x) <= 1
                                                else parser.error(f'CD-hit allows thresholds between 0.4 to 1'))
    parser.add_argument('--word_length', default='2', choices=['2', '3', '4', '5'],
                        help='A heuristic of CD-hit. Choose of word size:\n5 for similarity thresholds 0.7 ~ 1.0\n4 for similarity thresholds 0.6 ~ 0.7\n3 for similarity thresholds 0.5 ~ 0.6\n2 for similarity thresholds 0.4 ~ 0.5')
    parser.add_argument('--discard', default='1', help='Include only sequences longer than <$discard> for the analysis. (CD-hit uses only sequences that are longer than 10 amino acids. When the analysis includes shorter sequences, this threshold should be lowered. Thus, it is set here to 1 by default.)')
    parser.add_argument('--cluster_algorithm_mode', default='0', help='0 - clustered to the first cluster that meet the threshold (fast). 1 - clustered to the most similar cluster (slow)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
    
    start_time = time()
    cluster_sequences(args.fasta_file, args.output_prefix, args.done_file_path,
                      args.threshold, args.word_length, args.discard, args.cluster_algorithm_mode, sys.argv)
    end_time = time()
    logger.info(f'Done cluster sample name {args.sample_name} (took {measure_time(int(end_time-start_time))}).')
