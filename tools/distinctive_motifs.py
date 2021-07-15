
'''
Extract ranked distinctive motifs ignoring artifacts
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys 
import logging
from matplotlib.patches import Patch
import os
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)
from auxiliaries.pipeline_auxiliaries import log_scale


colors_map = {
    'perfect': 'green',
    'mixed': 'blue',
    'artifact': 'red',
    'negative': 'orange'
}

def get_sorted_features(feature_importance_path: str):
    with open(feature_importance_path, 'r') as f:
        features = [line.strip().split('\t') for line in f.readlines()]
    features = sorted(features, key=lambda x: float(x[1]), reverse=True)
    return features

def is_artifact(motif: str, df: pd.DataFrame, bio_cond: str, invalid_mix: str, score: float, order: int, logger: logging.Logger, txt_file_out: str):
    data = df[['label', 'sample_name', motif]]
    other_max = data.loc[data['label'] == 'other', motif].max()
    bc_min = data.loc[data['label'] == bio_cond, motif].min()
    artifact = bc_min < other_max
    is_perfect = bc_min > other_max
    is_valid_mix = True
    mixed_samples = list(data.loc[(data['label'] == 'other') & (data[motif] >= bc_min), 'sample_name'])
    if artifact:
        message = f'Motif {motif} is artifact, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n'
        txt_file_out.write(message) if txt_file_out else logger.info(message)
    elif is_perfect:
        message = f'Motif {motif} is perfect, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n'
        txt_file_out.write(message) if txt_file_out else logger.info(message)
    else:
        if invalid_mix:
            is_valid_mix = not any(invalid_mix in s for s in mixed_samples)
        if is_valid_mix:
            message = f'Motif {motif} is mixed, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n'
            txt_file_out.write(message) if txt_file_out else logger.info(message)
        else:
            message = f'Motif {motif} is invalid mix, other_max={other_max}, bc_min={bc_min}, importance score={score}, importance order={order}\n'
            txt_file_out.write(message) if txt_file_out else logger.info(message)
        # TODO convert mixed samples to mixed bc inc. count
    return artifact, is_perfect, is_valid_mix, mixed_samples


def generate_heatmap(base_path: str, df: pd.DataFrame, colors, title: str , rank_method: str, logger: logging.Logger):
    logger.info('Generating heatmap...')
    df.set_index('sample_name', inplace=True)
    df = log_scale(df, rank_method)

    map_path = f'{base_path}.svg'
    number_of_samples = df.shape[0]

    map = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True, col_colors=colors)
    plt.setp(map.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150 / number_of_samples)
    handles = [Patch(facecolor=colors_map[name]) for name in colors_map]
    plt.legend(handles, colors_map, title='Type', bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    map.ax_heatmap.set_title(title, pad=25, fontsize=14)
    map.savefig(map_path, format='svg', bbox_inches="tight")
    plt.close()


def save_output(base_path: str, data, logger: logging.Logger):
    logger.info('Saving results...')
    output_path = f'{base_path}.csv'
    df = pd.DataFrame(data, columns=['motif', 'label', 'is_artifact', 'is_perfect', 'is_valid_mix', 'mixed_samples', 'importance', 'order'])
    df.to_csv(output_path, index=False)


def extract_distinctive_motifs(data_path: str, count: int, feature_importance_path: str, biological_condition: str, output_base_path: str, title_heatmap: str, 
                           invalid_mix: str, epsilon: float, min_importance_score: float, rank_method: str, logger: logging.Logger):
    logger.info("Loading features...")
    features = get_sorted_features(feature_importance_path)
    data_df = pd.read_csv(data_path)
    
    total = 0
    last_score = 0
    last_order = 0
    distinctive_motifs = []
    i = 0
    perfect_count = 0
    artifact_count = 0
    invalid_mix_count = 0
    mixed_count = {}
    colors = []
    output = []
    output_label = ''
    txt_file_out = open(f'{output_base_path}.txt', 'w') if output_base_path else None
    for feature in features:
        i += 1
        motif = feature[0]
        score = float(feature[1])
        if score < min_importance_score:
            break
        artifact, is_perfect, is_valid_mix, mixed_samples = is_artifact(motif, data_df, biological_condition, invalid_mix, score, i, logger, txt_file_out)
        if total >= count and score + epsilon < last_score:
            break
        if artifact:
            artifact_count += 1
            colors.append(colors_map['artifact'])
            output_label = 'artifact'
            output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
            continue
        if not is_valid_mix:
            invalid_mix_count += 1
            colors.append(colors_map['negative'])
            output_label = 'negative'
            output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
            continue
        if is_perfect:
            perfect_count += 1
            colors.append(colors_map['perfect'])
            output_label = 'perfect'
        else:
            colors.append(colors_map['mixed'])
            output_label = 'mixed'
            for sample in mixed_samples:
                sample_count = 0
                try:
                    sample_count = mixed_count[sample]
                except:
                    pass
                mixed_count[sample] = sample_count + 1
        distinctive_motifs.append(motif)
        output.append([motif, output_label, artifact, is_perfect, is_valid_mix, mixed_samples, score, i])
        last_order = i
        total += 1
        if total == count:
            last_score = score
    
    summary = f'\nDistinctive motifs ({len(distinctive_motifs)}/{last_order} tested): {distinctive_motifs}\n' + \
              f'Perfects count: {perfect_count}\n' + \
              f'Mixed count: {mixed_count}\n' + \
              f'Filtered: Artifacts count: {artifact_count}, Invalid mixes count: {invalid_mix_count}\n'

    if output_base_path:
        txt_file_out.write(summary)
        txt_file_out.close()

        columns = ['sample_name'] + [x[0] for x in features[:count]]
        generate_heatmap(output_base_path, data_df[columns], colors, title_heatmap, rank_method, logger)
        save_output(output_base_path, output, logger)
    else:
        logger.info(summary)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A path for motif ranking CSV')
    parser.add_argument('count', type=int, help='Number of motifs to extract, e.g. 10 for top 10') 
    parser.add_argument('feature_importance_path', type=str, help='A path for the feature importance path for specific biological condition')
    parser.add_argument('biological_condition', type=str, help='Biological condition name')
    parser.add_argument('title_heatmap', type=str, help='A title for generated heatmap')
    parser.add_argument('--output_base_path', type=str, help='An optional base path for generated files, extensions will be added to the base path')
    parser.add_argument('--invalid_mix',type=str, default=None, help='A argument to know if there is compare to naive')
    parser.add_argument('--epsilon', type=int, default=0, help='range of mistake')
    parser.add_argument('--min_importance_score',type=int, default=0, help='The minimun score of motifs to take')
    parser.add_argument('--rank_method', type=str, choices=['pval', 'tfidf', 'shuffles', 'hits'], default='hits', help='Motifs ranking method')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    extract_distinctive_motifs(args.data_path, args.count, args.feature_importance_path, args.biological_condition, args.output_base_path, args.title_heatmap,
                                args.invalid_mix, args.epsilon, args.min_importance_score, args.rank_method, logger)
    