import logging
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, call



def aggregate_scores(scores_path, bc):

    # scores_path is a folder in which each file contains the scores of one of the scans split, e.g.:
    # /groups/pupko/orenavr2/igomeProfilingPipeline/experiments/test/analysis/model_fitting/17b/hits_scores
    output_path = f'{os.path.split(scores_path)[0]}/hits.txt'

    call(f'cat {scores_path}/*{bc}_motifs_*.txt > {output_path}', shell=True)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_path)


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('scores_path', type=lambda x: x.rstrip('/'),
                        help='A path to a (MEME) file with the motifs against which the random peptides were scanned (in silico)')
    parser.add_argument('biological_condition', help='The bc that should be aggregated')
    parser.add_argument('--new_rf', action='store_true', help='run new random forest version')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    aggregate_scores(args.scores_path, args.biological_condition)
