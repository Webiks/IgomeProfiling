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
from time import time


def calculate_pssm_thresholds(meme_path, cutoffs_path, faa_path, number_of_random_pssms, 
                              shuffles, shuffles_percent, shuffles_digits, rpm_factor, output_path, done_path, rank_method, argv='no_argv',
                              pssm_score_peptide='./PSSM_score_Peptide/PSSM_score_Peptide'):

    if not os.path.exists(output_path):
        # TODO: any modules to load?
        if rank_method == 'pval':
            cmd = f'{pssm_score_peptide} -pssm {meme_path} -pssm_cutoffs {cutoffs_path} -seq {faa_path} ' \
                f'-out {output_path} -NrandPSSM {number_of_random_pssms} -CalcPSSM_Pval'
            if rpm_factor:
                cmd += ' -useFactor'
            logger.info(f'{datetime.datetime.now()}: starting CalcPSSM_Pval. Executed command is:\n{cmd}')
        elif rank_method == 'tfidf':
            cmd = f'./hits_cpp/hits -m {meme_path} -c {cutoffs_path} -s {faa_path} -o {output_path} --outputSequences'
            logger.info(f'{datetime.datetime.now()}: starting TF-IDF\' hits. Executed command is:\n{cmd}')
        else:  # shuffles
            cmd = f'./hits_cpp/hits -m {meme_path} -c {cutoffs_path} -s {faa_path} -o {output_path} --shuffles {shuffles} '\
                f'--shufflesPercent {shuffles_percent} --shufflesDigits {shuffles_digits}'
            logger.info(f'{datetime.datetime.now()}: starting Shuffles\' hits. Executed command is:\n{cmd}')
            if rpm_factor:
                cmd += ' --useFactor'
            logger.info(f'{datetime.datetime.now()}: starting Shuffles\' hits. Executed command is:\n{cmd}')

        subprocess.run(cmd, shell=True)
    else:
        logger.info(f'{datetime.datetime.now()}: skipping scanning calculation as it is already exist at:\n{output_path}')

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file_path', help='A path to a meme file with motifs against which a set of random peptides will be scanned')
    parser.add_argument('cutoffs_file_path', help='A path to a cutoffs file (peptide above cutoff? -> peptide is part of the motif')
    parser.add_argument('faa_file_path', help='A path to a faa file with peptides to scan against the pssms in the meme file')
    parser.add_argument('rank_method', choices=['pval', 'tfidf', 'shuffles'], default='pval', help='Motifs ranking method')
    parser.add_argument('number_of_random_pssms', type=int, help='Number of pssm permutations')
    parser.add_argument('output_path', help='A path to which the Pvalues will be written to')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('--shuffles', default=5, type=int, help='Number of controlled shuffles permutations')
    parser.add_argument('--shuffles_percent', default=0.2, type=float, help='Percent from shuffle with greatest number of hits (0-1)')
    parser.add_argument('--shuffles_digits', default=2, type=int, help='Number of digits after the point to print in scanning files')
    parser.add_argument('--no_rpm_factor', action='store_false', help='Disable multiplication hits by factor rpm for normalization')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    start = time()
    calculate_pssm_thresholds(args.meme_file_path, args.cutoffs_file_path, args.faa_file_path,
                              args.number_of_random_pssms, args.shuffles, args.shuffles_percent, args.shuffles_digits, args.no_rpm_factor,
                              args.output_path, args.done_file_path, args.rank_method, argv=sys.argv)
    end = time()
    print(f'total time (sec): {end - start}')
