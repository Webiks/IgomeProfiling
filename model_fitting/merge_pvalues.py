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


def get_results_pval(consensusesValues, sample_name, path):
    # pvalues = []
    # hits = []
    with open(path) as f:
        for line in f:
            if line.startswith('##'):
                # "## PSSM_name	p_Value	True_Hits: num_of_hits"
                continue
            # "CLKGASFLAC_17b_clusterRank_0_uniqueMembers_339_clusterSize_2659173.71.faa\t0.01\tTrue_Hits: 118"
            motif = line.split('_')[0]
            line_tokens = line.split('\t')
            pvalue = line_tokens[1]
            hits = line_tokens[2].split()[-1]

            if motif not in consensusesValues:
                motifSamples = {}
                consensusesValues[motif] = motifSamples
            else:
                motifSamples = consensusesValues[motif]

            motifSamples[sample_name] = { 'hits': hits, 'pvalue': pvalue }


def get_results_shuffles(consensusesValues, sample_name, path):
    with open(path) as f:
        for line in f:
            # first line get motif - MOTIF CDWFEQYGLRLR_17b_clusterRank_0003_uniqueMembers_51_clusterSize_62065.05.faa
            motif = line.split()[1].split('_')[0]
            # secand line get hits -  HITS 1043
            hits = f.readline().split()[1]
            #third line skip
            line3 = f.readline()
            # four line get valeus - RANK 1.00
            pvalue = f.readline().split()[1]

            if motif not in consensusesValues:
                motifSamples = {}
                consensusesValues[motif] = motifSamples
            else:
                motifSamples = consensusesValues[motif]

            motifSamples[sample_name] = { 'hits': hits, 'pvalue': pvalue }


def aggregate_values_results(meme_path, scanning_results_dir_path, bc, samplename2biologicalcondition_path,
                              aggregated_pvalues_path, aggregated_hits_path, done_path, rank_method, bc_sample_names, argv='no_argv'):
    if not bc_sample_names:
        samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path,
                                                'Barcode {} belongs to more than one sample_name!!')
        bc_sample_names = [sample for sample in samplename2biologicalcondition if samplename2biologicalcondition[sample] == bc]
    
    all_consensuses = get_consensus_sequences_from_meme(meme_path)
    samples = set()
    # consensuses => samples => hits, pvals
    consensusesValues = {}
    for file_name in sorted(os.listdir(scanning_results_dir_path)):
        sample_name = file_name.split('_peptides')[0]
        samples.add(sample_name)
        if rank_method == 'pval':
            get_results_pval(consensusesValues, sample_name, os.path.join(scanning_results_dir_path, file_name))
        else: # rank_method == 'shuffles'
            get_results_shuffles(consensusesValues, sample_name, os.path.join(scanning_results_dir_path, file_name))
    pvalues_f = open(aggregated_pvalues_path, 'w')
    hits_f = open(aggregated_hits_path, 'w')

    # header
    header = f'sample_name,label,{",".join(all_consensuses).rstrip(",")}\n'
    pvalues_f.write(header)
    hits_f.write(header)

    for sample in sorted(samples):
        label = bc if sample in bc_sample_names else 'other'
        pvalues_f.write(f'{sample},{label}')
        hits_f.write(f'{sample},{label}')
        for consensus in all_consensuses:
            values = consensusesValues[consensus][sample]
            pvalues_f.write(f',{values["pvalue"]}')
            hits_f.write(f',{values["hits"]}')
        pvalues_f.write('\n')
        hits_f.write('\n')

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
    parser.add_argument('biological_condition', help='Positive class\' label. All samples of another biological condition will be labeled as "other"')
    parser.add_argument('aggregated_pvalues_path', help='A path to which the Pvalues table will be written to')
    parser.add_argument('aggregated_hits_path', help='A path to which the hits table will be written to')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('--rank_method', choices=['pval', 'shuffles'], default='shuffles', help='Motifs ranking method')
    parser.add_argument('--bc_sample_names', type=str, help='Samples that are for the label bc, seperate by comma')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    aggregate_values_results(args.meme_path, args.scanning_results_dir_path,
                            args.biological_condition, args.samplename2biologicalcondition_path,
                            args.aggregated_pvalues_path, args.aggregated_hits_path,
                            args.done_file_path, args.rank_method, args.bc_sample_names, sys.argv)

