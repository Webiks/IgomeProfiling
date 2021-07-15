import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, remove_redundant_newlines_from_fasta, \
    get_unique_members_from


def reconstruct_msa(sequences_file_path, output_file_path, done_path, argv='no_argv'):
    number_of_unique_members = get_unique_members_from(sequences_file_path)
    if number_of_unique_members > 1:
        import subprocess
        # TODO: module load mafft..
        # --auto Automatically selects an appropriate strategy from L-INS-i, FFT-NS-i and FFT-NS-2, according to data size.
        # --amino tells mafft that's an amino acid msa. If you let it decide by itself, it might wrong on small data sets
        # as they might look like dna but they are NOT! e.g.,
        # [orenavr2@powerlogin-be2 test]$ cat /groups/pupko/orenavr2/igomeProfilingPipeline/experiments/test/analysis/motif_inference/17b_03/unaligned_sequences/17b_03_clusterRank_215_uniqueMembers_2_clusterSize_252.81.faa
        # >seq_235_lib_12_len_12_counts_126.40626975097965
        # CNTDVACAAPGN
        # >seq_1112_lib_C8C_len_10_counts_126.40626975097965
        # CTTACAPVNC
        cmd = f'mafft --auto --amino {sequences_file_path} > {output_file_path}'
        logger.info(f'{datetime.datetime.now()}: Starting MAFFT. Executed command is:\n{cmd}')
        subprocess.run(cmd, shell=True)
    else:
        logger.info(f'{datetime.datetime.now()}: skipping alignment for a cluster with a single member. '
                    f'Writing the output file as is to\n'
                    f'{output_file_path}')
        with open(sequences_file_path) as unaligned_f:
            content = unaligned_f.read()
        with open(output_file_path, 'w') as aligned_f:
            aligned_f.write(content)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_file_path)

    # override the results with clean ones (no redundant new lines. For further details see function's doc)
    remove_redundant_newlines_from_fasta(output_file_path, output_file_path)

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(output_file_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sequences_file_path', help='A path to a file with unaligned sequences')
    parser.add_argument('output_file_path', help='A path to a file in which the aligned sequences will be written')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    reconstruct_msa(args.sequences_file_path, args.output_file_path, args.done_file_path, sys.argv)
