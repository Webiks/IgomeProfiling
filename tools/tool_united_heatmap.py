import sys
import pandas as pd
import os
import logging
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import load_table_to_dict, log_scale


def get_motifs_importance(biological_condition, bc_dir_path, rank_type, min_num_motifs, max_num_motifs, max_difference_from_last_motif, max_difference_from_fitst_motif):
    important_motifs = []
    motif_importance_path = os.path.join(bc_dir_path, f'{biological_condition}_{rank_type}_model', 'best_model', 'feature_importance.txt')
    dict_importance = load_table_to_dict(motif_importance_path, 'Motif {} is not unique!!')
    first_val = None
    last_val = 0
    for key, value in dict_importance.items():
        if first_val is None:
            first_val = float(value)
            last_val = first_val
            important_motifs.append(key)
            continue 
        current_val = float(value)
        if len(important_motifs) >= max_num_motifs or first_val - current_val >=  max_difference_from_fitst_motif or last_val - current_val >= max_difference_from_last_motif:
            break
        important_motifs.append(key)
        last_val = current_val
    if len(important_motifs) <= min_num_motifs and len(dict_importance.keys())>= min_num_motifs:
        important_motifs = list(dict_importance.keys())[:min_num_motifs]
    return important_motifs

def united_csv(input_path, samplename2biologicalcondition_path, min_num_motifs, max_num_motifs, max_difference_from_last_motif, max_difference_from_fitst_motif):
    samplename2biologicalcondition = load_table_to_dict(samplename2biologicalcondition_path, 'Barcode {} belongs to more than one sample_name!!')
    biological_conditions = sorted(set(samplename2biologicalcondition.values()))
    color_name = sns.color_palette(None, len(biological_conditions))
    dict_color = dict(zip(biological_conditions, color_name))
    list_hits = []
    list_values = []
    color_motifs_hits = []
    color_motifs_values = []
    for bc, color in zip(biological_conditions, color_name):
        bc_dir_path = os.path.join(input_path, bc)
        list_motifs_hits = get_motifs_importance(bc, bc_dir_path, 'hits', min_num_motifs, max_num_motifs, max_difference_from_last_motif, max_difference_from_fitst_motif)
        list_motifs_values = get_motifs_importance(bc, bc_dir_path, 'values', min_num_motifs, max_num_motifs, max_difference_from_last_motif, max_difference_from_fitst_motif)
        bc_hits = os.path.join(bc_dir_path, f'{bc}_hits.csv')
        bc_values = os.path.join(bc_dir_path, f'{bc}_values.csv')
        df_hits = pd.read_csv(bc_hits, index_col=0)
        df_value = pd.read_csv(bc_values, index_col=0)
        df_hits = df_hits[list_motifs_hits]
        df_value = df_value[list_motifs_values]
        
        number_of_motifs_hits = len(list_motifs_hits)
        number_of_motifs_values = len(list_motifs_values)
        color_motifs_hits.append([color] * number_of_motifs_hits)
        color_motifs_values.append([color] * number_of_motifs_values)
        list_hits.append(df_hits)
        list_values.append(df_value)
    #merge all the df to one df. 
    hits_all=pd.concat(list_hits, axis=1) 
    values_all = pd.concat(list_values, axis=1)
    hits_all = hits_all.reset_index()
    values_all = values_all.reset_index()

    color_motifs_hits = sum(color_motifs_hits, [])
    color_motifs_values = sum(color_motifs_values, [])
    return hits_all, values_all, color_motifs_hits, color_motifs_values, dict_color

def generate_heat_map(df, rank_method, number_of_samples, output_path, color_list, dict_color):
    df = df.set_index(df.columns[0])
    train_data = log_scale(df, rank_method)
    cm = sns.clustermap(train_data, cmap="Blues", col_cluster=False, yticklabels=True, col_colors=color_list)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    handles = [Patch(facecolor=dict_color[bc]) for bc in dict_color]
    plt.legend(handles, dict_color, title='BC',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    cm.ax_heatmap.set_title(f"A heat-map of the all biological condition discriminatory motifs", pad=25.0)
    cm.savefig(f"{output_path}.svg", format='svg', bbox_inches="tight")
    plt.close()

def united_heatmap(data_path, output_path, samplename2biologicalcondition_path, min_num_motifs, max_num_motifs, max_difference_from_last_motif, max_difference_from_fitst_motif, rank_method):
    
    hits, values, color_motifs_hits, color_motifs_values, dict_color = united_csv(data_path, samplename2biologicalcondition_path, min_num_motifs, max_num_motifs,
                                                                                  max_difference_from_last_motif, max_difference_from_fitst_motif)
    output_path_hits = os.path.join(output_path,'hits_all_bc')
    output_path_values = os.path.join(output_path,'values_all_bc')
    generate_heat_map(hits, 'hits', len(hits), output_path_hits, color_motifs_hits, dict_color)
    generate_heat_map(values, rank_method, len(values), output_path_values, color_motifs_values, dict_color)

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A path for all the csv files with data matrix to model ')
    parser.add_argument('output_path', type=str, help='A path for folder to put the output')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='Path for a file samplename2biologicalcondition')
    parser.add_argument('--min_num_motifs', type=int, default=1, help='Minimum number of motifs to united from every BC')
    parser.add_argument('--max_num_motifs', type=int, default=10, help='Maximum number of motifs to united from every BC')
    parser.add_argument('--max_difference_from_last_motif', type=float, default=0.01, help='Take the motif if the difference of his importent values is less than the last motif')
    parser.add_argument('--max_difference_from_fitst_motif', type=float, default=0.05, help='Take the motif if the difference of his importent values is less than the first motif')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles'], default='shuffles', help='Rank method')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
    
    united_heatmap(args.data_path, args.output_path, args.samplename2biologicalcondition_path, args.min_num_motifs, args.max_num_motifs, 
                    args.max_difference_from_last_motif, args.max_difference_from_fitst_motif, args.rank_method)
