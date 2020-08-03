import logging
import os
import sys
import pandas as pd

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


def get_results(consensusesValues, sample_name, path):
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


def aggregate_pvalues_results(meme_path, scanning_results_dir_path, bc, samplename2biologicalcondition_path,
                              aggregated_pvalues_path, aggregated_hits_path, use_new_rf, done_path, argv='no_argv'):

    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path,
                                                'Barcode {} belongs to more than one sample_name!!')

    all_consensuses = get_consensus_sequences_from_meme(meme_path)

    if use_new_rf:
        pvalues_f = open(aggregated_pvalues_path, 'w')
        hits_f = open(aggregated_hits_path, 'w')

        #header
        pvalues_result = hits_result = f'sample_name,label,{",".join(all_consensuses)}'
        for file_name in sorted(os.listdir(scanning_results_dir_path)):
            if file_name.endswith('100.txt'):
                raise TypeError  # why?
            if file_name.startswith('.'):
                # system file...
                continue

            if file_name.endswith('00.txt'):
                # next sample is starting
                pvalues_f.write(f'{pvalues_result.rstrip(",")}\n')
                hits_f.write(f'{hits_result.rstrip(",")}\n')
                sample_name = file_name.split('_peptides')[0]
                if bc in sample_name:
                    label = samplename2biologicalcondition[sample_name]
                else:
                    label = 'other'
                pvalues_result = hits_result = f'{sample_name},{label},'

            pvalues, hits = get_results(os.path.join(scanning_results_dir_path, file_name))
            pvalues_result += ','.join(pvalues) + ','
            hits_result += ','.join(hits) + ','

        pvalues_f.write(f'{pvalues_result.rstrip(",")}\n')
        hits_f.write(f'{hits_result.rstrip(",")}\n')

        pvalues_f.close()
        hits_f.close()

        # remove insignificant features:
        df = pd.read_csv(aggregated_pvalues_path)
        # features with at least one significant score, across positive-labeled samples
        positive_class_df = df[df['label'] != 'other']
        significant_features = (positive_class_df.drop(['sample_name', 'label'], axis=1) < 0.05).sum()>0
        mask = pd.concat([pd.Series([True, True], index=['sample_name', 'label']),
                        significant_features])
        df = df.loc[:, mask]
        # df = pd.concat([df.loc[:, ['sample_name', 'label']], df.drop(['sample_name', 'label'], axis=1).loc[:, significant_features]], axis=1)
        df.to_csv(aggregated_pvalues_path.replace('_insignificant', ''), index=False)
    
    else:
        samples = set()
        # consensuses => samples => hits, pvals
        consensusesValues = {}
        for file_name in sorted(os.listdir(scanning_results_dir_path)):
            sample_name = file_name.split('_peptides')[0]
            samples.add(sample_name)
            get_results(consensusesValues, sample_name, os.path.join(scanning_results_dir_path, file_name))

        pvalues_f = open(aggregated_pvalues_path, 'w')
        hits_f = open(aggregated_hits_path, 'w')

        # header
        header = f'sample_name,label,{",".join(all_consensuses).rstrip(",")}\n'
        pvalues_f.write(header)
        hits_f.write(header)

        for sample in sorted(samples):
            label = bc if samplename2biologicalcondition[sample] == bc else 'other'
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
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('--new_rf', action='store_true', help='run new random forest version')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    aggregate_pvalues_results(args.meme_path, args.scanning_results_dir_path,
                              args.biological_condition, args.samplename2biologicalcondition_path,
                              args.aggregated_pvalues_path, args.aggregated_hits_path, args.new_rf,
                              args.done_file_path, argv=sys.argv)


