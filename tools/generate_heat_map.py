'''
Extract ranked distinctive motifs ignoring artifacts
'''
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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


def process(values_path: str, motifs, colors, heatmap_title: str, output_base_path: str):
    values = pd.read_csv(values_path)
    columns = ['sample_name'] + motifs
    generate_heatmap(output_base_path, values[columns], colors, heatmap_title)


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    values_path = path.join(base_path, 'Texas_2012/Texas_2012_values_exp4_complete.csv')
    output_base_path = path.join(base_path, 'Texas_2012/cross_final_exp4_complete')
    heatmap_title = 'Exp8 cross best motifs on Exp4'
    # Green - top10 and perfect
    # Blue - top10 and mixed,
    # LightGreen - unexpected perfect
    # Orange - unexpected and mixed
    motifs = ['ADLAFRWGGM', 'SLFLHELTKSSG', 'LYFTNPPEPC', 'FQEWDMTRHS', 'WLADPSRVRGT', 'SRSAFLYPPP', 'CMNFNQNSSC', 'ARSPTNANIS', 'CLEWPSSIANC']
    colors = ['green', 'green', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange']
    process(values_path, motifs, colors, heatmap_title, output_base_path)

