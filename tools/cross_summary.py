"""
Union distinctive and identified motifs to single table
"""
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re


colors_map = {
    'perfect': 'green',
    'mixed': 'blue',
    'incorrect': 'red',
    'unexpected': 'orange',
    'artifact': 'red',
    'negative': 'orange',
    'irrelevant': 'grey'
}


cluster_pattern = re.compile('(.+?)_(.+)_clusterRank')


def load_union(path: str, base_name: str):
    motifs = {}
    with open(path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        clusters = line.split(',')
        clustered = {}
        for cluster in clusters:
            cluster_match = cluster_pattern.match(cluster)
            consensus = cluster_match.group(1)
            biological_condition = cluster_match.group(2)
            clustered[consensus] = 'base' if biological_condition == base_name else 'other'
        for key in clustered:
            union = { 'base': [], 'other': [] }
            for other_key in clustered:
                if key == other_key: continue
                union[clustered[other_key]].append(other_key)
            motifs[key] = union
    return motifs


def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0, color_from_index=2,
                     color_to_index = 2, output_path=None, ax=None, has_union_data=None,
                     **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    cell_align = 'center' if has_union_data else 'left'
    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, cellLoc=cell_align, **kwargs)
    if has_union_data:
        mpl_table.auto_set_column_width(col=list(range(len(data.columns))))
    else:
        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)
    
    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        elif k[1] >= color_from_index and k[1] < color_to_index:
            cell.set_facecolor(colors_map[data.iat[k[0] - 1, k[1]]])
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    if output_path:
        fig.savefig(output_path)
    return ax.get_figure(), ax

# TODO refactor
def process(results, top_motifs, union_data, output_base_path: str):
    motifs = {}
    # Index
    for result in results:
        source = result[0]
        df = pd.read_csv(result[1])
        df = df[['motif', 'label']]
        for _, row in df.iterrows():
            motif = row['motif']
            label = row['label']          
            motif_data = None
            try:
                motif_data = motifs[motif]
            except:
                motif_data = { 'is_top': motif in top_motifs }
                motifs[motif] = motif_data
            motif_data[source] = label
    # TODO extract details (mixes)? maybe with flag
    # Flatten
    motif_keys = list(motifs.keys())
    first_key = motif_keys[0]
    columns = list(motifs[first_key].keys())
    data = []
    for motif in motif_keys:
        motif_data = motifs[motif]
        motif_flat_data = [motif]
        for column in columns: motif_flat_data.append(motif_data[column])
        if union_data:
            motif_flat_data.append(', '.join(union_data[motif]['base']))
            motif_flat_data.append(', '.join(union_data[motif]['other']))
        data.append(motif_flat_data)
    # Save to csv
    print('Saving results...')
    output_path = f'{output_base_path}.csv'
    cols = ['motif'] + columns
    if union_data: cols += ['self_union', 'other_union']
    df = pd.DataFrame(data, columns=cols)
    df.to_csv(output_path, index=False)
    # Save table (svg)
    output_path = f'{output_base_path}_all.svg'
    has_union_data = union_data is not None
    render_mpl_table(df, output_path=output_path, color_to_index=len(columns) + 1, has_union_data=has_union_data)
    df_groupped = df.groupby('is_top')
    for is_top, group in df_groupped:
        name = 'top' if is_top else 'others'
        output_path = f'{output_base_path}_{name}.svg'
        group = group.drop('is_top', axis=1)
        render_mpl_table(group, output_path=output_path, color_from_index=1, color_to_index=len(columns), has_union_data=has_union_data)


if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting'
    top_motifs = ['CSTTTTTRTC', 'CDTDGWIRMDPC', 'CPSHISLRNPC', 'CTRPTRPSLWMSPAC', 'CTNLCGTLAC', 'CCHHIFQRPPC', 'CPACSLFTC', 'CDAHSLFFPC', 'CGGTAGSFSRC', 'CLSSLFTRC']
    results = [
        ['base', path.join(base_path, 'ferret_texas/distinctive_motifs_all.csv')],
        ['reduced', path.join(base_path,'ferret_texas/exp4_texas_motifs_on_exp8_reduced.csv')],
        ['complete', path.join(base_path,'ferret_texas/exp4_texas_motifs_on_exp8_complete.csv')],
    ]
    base_name = 'ferret_texas'
    union_data = load_union('/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_inference/cross/texas_combined.csv', base_name)
    output_base_path = path.join(base_path, 'ferret_texas/cross_summary')
    process(results, top_motifs, None, output_base_path)
    process(results, top_motifs, union_data, output_base_path + "_union")