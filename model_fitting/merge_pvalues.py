import logging
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
else:
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import verify_file_is_not_empty, load_table_to_dict


def get_consensus_sequences_from_meme(meme_path):
    result = []
    with open(meme_path) as f:
        for line in f:
            if not line.startswith('MOTIF'):
                continue
            # "MOTIF CLKGASFLAC_17b_clusterRank_0_uniqueMembers_339_clusterSize_2659173.71.faa"
            result.append(line.split()[1].split('_')[0])

    return result


def get_results(path):
    pvalues = []
    hits = []
    with open(path) as f:
        for line in f:
            if line.startswith('##'):
                # "## PSSM_name	p_Value	True_Hits: num_of_hits"
                continue
            # "CLKGASFLAC_17b_clusterRank_0_uniqueMembers_339_clusterSize_2659173.71.faa\t0.01\tTrue_Hits: 118"
            line_tokens = line.split('\t')
            # consensus.append(line_tokens[0].split('_')[0])
            pvalues.append(line_tokens[1])
            hits.append(line_tokens[2].split()[-1])

    return pvalues, hits


def aggregate_pvalues_results(meme_path, scanning_results_dir_path, samplename2biologicalcondition_path,
                              aggregated_pvalues_path, aggregated_hits_path, done_path, argv='no_argv'):

    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path,
                                                'Barcode {} belongs to more than one sample_name!!')

    all_consensuses = get_consensus_sequences_from_meme(meme_path)
    pvalues_f = open(aggregated_pvalues_path, 'w')
    hits_f = open(aggregated_hits_path, 'w')

    #header
    pvalues_result = hits_result = f'sample_name,label,{",".join(all_consensuses)}'
    # consensus2pvalues_column = {consensus:[] for consensus in all_consensuses}
    # consensus2hits_column = {consensus:[] for consensus in all_consensuses}
    for file_name in sorted(os.listdir(scanning_results_dir_path)):
        if file_name.endswith('100.txt'):
            raise TypeError

        if file_name.endswith('00.txt'):
            # next sample is starting
            pvalues_f.write(f'{pvalues_result.rstrip(",")}\n')
            hits_f.write(f'{hits_result.rstrip(",")}\n')
            sample_name = file_name.split('_peptides')[0]
            label = samplename2biologicalcondition[sample_name]
            pvalues_result = hits_result = f'{sample_name},{label},'

        pvalues, hits = get_results(os.path.join(scanning_results_dir_path, file_name))
        pvalues_result += ','.join(pvalues) + ','
        hits_result += ','.join(hits) + ','

    pvalues_f.write(f'{pvalues_result.rstrip(",")}\n')
    hits_f.write(f'{hits_result.rstrip(",")}\n')

    pvalues_f.close()
    hits_f.close()

    # make sure that there are results and the file is not empty
    verify_file_is_not_empty(aggregated_pvalues_path)
    verify_file_is_not_empty(aggregated_hits_path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('meme_path', help='A path to a (MEME) file with the motifs against which the random peptides were scanned (in silico)')
    parser.add_argument('scanning_results_dir_path', help='A path in which each file contains a hits/pvals computation. '
                                                         'These files will be aggregated into one table.')
    parser.add_argument('aggregated_pvalues_path', help='A path to which the Pvalues table will be written to')
    parser.add_argument('aggregated_hits_path', help='A path to which the hits table will be written to')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    aggregate_pvalues_results(args.meme_path, args.scanning_results_dir_path,
                              args.samplename2biologicalcondition_path,
                              args.aggregated_pvalues_path, args.aggregated_hits_path,
                              args.done_file_path, argv=sys.argv)

