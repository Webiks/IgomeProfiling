import datetime
import os
import sys
from collections import defaultdict

if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)
from auxiliaries.pipeline_auxiliaries import load_table_to_dict

import logging
logger = logging.getLogger('main')


def merge_meme_files(motif_inference_path, biological_condition, merged_meme_path, done_path, samplename2biologicalcondition_path, samples_to_skip, argv='no_argv'):
    """
    :param motif_inference_path: A path in which each folder corresponds to a sample and contains a meme file for the
    motifs in this sample.
    :param biological_condition: A biological condition to merge its samples meme files
    :param meme_file_name: An optional file name for the output file
    :param samplename2biologicalcondition_path: A file that connect a sample name to is biological condition
    :return: A merged meme file at $motif_inference_path/$biological_condition/$meme_file_name
    """
    logger.info(f'{datetime.datetime.now()}: merging meme files of {biological_condition}')

    memes = []
    first_meme = True
    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path,
                                                'Barcode {} belongs to more than one sample_name!!')
    list_sample_of_bc=[i for i in samplename2biologicalcondition if samplename2biologicalcondition[i]==biological_condition]
    for sample_name in sorted(list_sample_of_bc):  # sample name of the specific bc
        dir_path = os.path.join(motif_inference_path, sample_name)
        if not os.path.isdir(dir_path):
            # skip file or folders of non-related biological condition
            continue
        if sample_name in samples_to_skip:
            # skip unwanted samples
            logger.info(f'{datetime.datetime.now()}: skipping {dir_path} (according to user\'s request)')
            continue
        try:
            with open(os.path.join(dir_path, 'meme.txt')) as f:
                if first_meme:
                    meme_header = ''.join(f.readline() for i in range(6))  # header's length is 6 lines
                    first_meme = False
                else:
                    # skip meme header (6 rows)
                    for i in range(6):
                        f.readline()
                # add current memes to memes list
                memes.extend(f.read().rstrip().split('\n\n\n'))
        except:
            logger.error(f'{datetime.datetime.now()}: no meme file in {dir_path}\n'
                         f'This is its content:\n{os.listdir(dir_path)}')

    # sort memes accross (biological condition) by their cluster size
    memes.sort(key=lambda meme: float(meme.split('.faa')[0].split('clusterSize_')[-1]), reverse=True)

    with open(merged_meme_path, 'w') as f:
        f.write(meme_header)
        f.write('\n\n\n'.join(memes))
        f.write('\n\n\n')

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')



if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('motif_inference_path', help='A path in which each folder corresponds to a sample and contains a meme file for the motifs in this sample.')
    parser.add_argument('biological_condition', help='A biological condition to merge its samples meme files')
    parser.add_argument('merged_meme_path', help='A path to the output file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('samplename2biologicalcondition_path', help='A path to the sample name to biological condition file.')
    parser.add_argument('--skip_sample', default='a_weird_str_that_shouldnt_be_a_sample_name_by_any_chance',
                        help='A sample name that should be skipped, e.g., for testing purposes. More than one sample '
                             'name should be separated by commas but no spaces. '
                             'For example: 17b_05,17b_05_test,another_one')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    merge_meme_files(args.motif_inference_path, args.biological_condition, args.merged_meme_path,
                     args.done_file_path, args.samplename2biologicalcondition_path ,args.skip_sample.split(','), sys.argv)