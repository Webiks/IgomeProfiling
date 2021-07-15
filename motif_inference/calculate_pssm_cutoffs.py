import datetime
import subprocess
import logging
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty


def calculate_pssm_cutoffs(meme_path, output_path, done_path, total_memes, cutoff_random_peptitdes_percentile, argv='no_argv',
                           pssm_score_peptide='./PSSM_score_Peptide/PSSM_score_Peptide'):

    if not os.path.exists(output_path):
        # TODO: any modules to load?
        cmd = f'{pssm_score_peptide} -pssm {meme_path} -pssm_cutoffs {output_path} -CalcPSSM_Cutoff -total_memes {total_memes} -cutoff_random_peptitdes_percentile {cutoff_random_peptitdes_percentile}'
        logger.info(f'{datetime.datetime.now()}: starting PSSM_score_Peptide. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
    else:
        logger.info(f'{datetime.datetime.now()}: skipping cutoffs calculation as it is already exist at:\n{output_path}')

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file_path', help='A path to a meme file')
    parser.add_argument('output_path', help='A path in which a new subfolder with the united motifs will be written to')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('--total_memes', type=int, default=0, help='Total memes in biological condition. Used if input data is splitted')
    parser.add_argument('--cutoff_random_peptitdes_percentile', type=float, default=0.05, help='Calculate cutoff (hit threshold) from random peptides\' top percentile score')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    calculate_pssm_cutoffs(args.meme_file_path, args.output_path, args.done_file_path, args.total_memes, args.cutoff_random_peptitdes_percentile, argv=sys.argv)
