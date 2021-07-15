"""
Extract identified biological conditions per motif
"""
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


colors_map = {
    'perfect': 'green',
    'mixed': 'blue',
    'incorrect': 'red',
    'unexpected': 'orange',
    'artifact': 'red',
    'negative': 'orange',
    'irrelevant': 'grey'
}


def load_samples_map(samples_to_bio_conds_path: str):
    with open(samples_to_bio_conds_path, 'r') as f:
        return dict([line.strip().split('\t') for line in f.readlines()])


def save_output(base_path: str, data):
    print('Saving results...')
    output_path = f'{base_path}.csv'
    df = pd.DataFrame(data, columns=['motif', 'label', 'is_perfect', 'is_positive', 'is_negative', 'identifies'])
    df.to_csv(output_path, index=False)


def generate_heatmap(base_path: str, df: pd.DataFrame, colors, title: str):
    print('Generating heatmap...')

    df.set_index('sample_name', inplace=True)
    map_path = f'{base_path}.svg'
    number_of_samples = df.shape[0]

    map = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True, col_colors=colors)
    plt.setp(map.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150 / number_of_samples)
    map.ax_heatmap.set_title(title, pad=25, fontsize=14)
    map.savefig(map_path, format='svg', bbox_inches="tight")
    plt.close()


def identify(values_path: str, samples_to_bio_conds_path: str, assumed_motifs, positive_bio_cond: str, negative_bio_cond: str, output_base_path: str, heatmap_title: str):
    samples_to_bio_conds = load_samples_map(samples_to_bio_conds_path)
    values = pd.read_csv(values_path)
    values['label'] = [samples_to_bio_conds[sample] for sample in values['sample_name']]
    motifs = list(values.columns)[2:]
    
    total_positive = 0
    total_negative = 0
    total_others = 0
    total_perfect = 0
    total_per_bio_cond = {}
    assumption_in = []
    assumption_out = []
    assumption_others = []
    output = []
    colors = []
    motif_label = ''
    for motif in motifs:
        is_negative = False
        is_positive = False
        is_perfect = False
        data = values[['label', motif]].groupby("label")
        scored_data = 2 * data[motif].max() + (1 - data[motif].std()) # TODO add configurable weights
        scored_data.sort_values(0, False, True)
        identified_bio_conds = []
        best_value = None
        for label, value in scored_data.iteritems():
            if not identified_bio_conds:
                best_value = value
                identified_bio_conds.append(label)
            elif value == best_value:
                identified_bio_conds.append(label)
            else:
                break
        if positive_bio_cond:
            is_positive = any(positive_bio_cond.lower() in s.lower() for s in identified_bio_conds)
        if negative_bio_cond:
            is_negative = any(negative_bio_cond.lower() in s.lower() for s in identified_bio_conds)
        if is_positive:
            total_positive += 1
            if len(identified_bio_conds) == 1: is_perfect = True
        if is_negative: total_negative += 1
        if is_perfect: total_perfect += 1
        if not is_positive and not is_negative: total_others += 1
        for bio_cond in identified_bio_conds:
            count = 0
            try:
                count = total_per_bio_cond[bio_cond]
            except:
                pass
            total_per_bio_cond[bio_cond] = count + 1
        if assumed_motifs:
            is_assumed = motif in assumed_motifs
            if is_assumed and is_positive and not is_negative:
                assumption_in.append(motif)
                motif_label = 'perfect' if is_perfect else 'mixed'
            elif is_assumed:
                assumption_out.append(motif)
                motif_label = 'incorrect'
            elif is_positive and not is_negative:
                assumption_others.append(motif)
                motif_label = 'unexpected'
            else:
                motif_label = 'irrelevant'
        else:
            if is_perfect:
                motif_label = 'perfect'
            elif is_positive and not is_negative:
                motif_label = 'mixed'
            elif is_negative and not is_positive:
                motif_label = 'negative'
            else:
                motif_label = 'artifact'
        output.append([motif, motif_label, is_perfect, is_positive, is_negative, identified_bio_conds])
        colors.append(colors_map[motif_label])
        print(f'{motif}: positive={is_positive}, negative={is_negative}, perfect={is_perfect}, identifies={identified_bio_conds}')
    
    print(f'\nTotal motifs={len(motifs)}, positives={total_positive}, negative={total_negative}, others={total_others}, perfects={total_perfect}')
    print(f'Count per biological condition: {total_per_bio_cond}')
    if assumed_motifs:
        print(f'\nAssumed count motifs={len(assumed_motifs)}, in={len(assumption_in)}, out={len(assumption_out)}, not expected={len(assumption_others)}')
        if assumption_in: print(f'Correctly assumed: {assumption_in}')
        if assumption_out: print(f'Incorrectly assumed: {assumption_out}')
        if assumption_others: print(f'Unexpected positives: {assumption_others}')
    if output_base_path:
        columns = ['sample_name'] + motifs
        generate_heatmap(output_base_path, values[columns], colors, heatmap_title)
        save_output(output_base_path, output)


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    assumed_motifs = ['CVGFVCTGPV', 'CNSNLSRSENL', 'FYSRTQYKFCEVGC', 'QTPYGRHPTL', 'YPKQDQPDIA', 'FFGALRP', 'CGSNYPRLSNI', 'CSGAAFEKYM', 'RDRPYMVDSPVN', 'YTRSSGMLIDC']
    values_path = path.join(base_path, 'Naive/Naive_values_exp4_complete.csv')
    hits_path = path.join(base_path, 'Naive/Naive_hits_exp4_complete.csv') # TODO check hits "back-up"
    samples_to_bio_conds_path = path.join(base_path, 'exp4+9_all_samples2bioconds.txt')
    positive_bio_cond = 'naive' # 'texas'
    negative_bio_cond = None # 'naive'
    output_base_path = path.join(base_path, 'Naive/exp8_naive_motifs_on_exp4_complete')
    heatmap_title = 'Exp8 Naive Motifs on Exp4 all samples'

    identify(values_path, samples_to_bio_conds_path, assumed_motifs, positive_bio_cond, negative_bio_cond, output_base_path, heatmap_title)