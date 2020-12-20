'''
Generate a heatmap
can call the function from the pipeline, or generate a specific heatmap and call the python file. 
'''
from os import path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

def generate_heatmap(output_path: str, df: pd.DataFrame, title: str, colors=False):
    print('Generating heatmap...')
    df=df.drop('label', 1) if 'label' in list(df.columns) else df
    df.set_index('sample_name', inplace=True)
    is_hits=True if 'hits' in output_path else False
    df=np.log2(df+1) if is_hits else df 
    
    map_path = f'{output_path}.svg'
    number_of_samples = df.shape[0]
    
    if colors: 
        map = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True,col_colors=colors)
    else:
        map = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(map.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150 / number_of_samples)
    map.ax_heatmap.set_title(title, pad=25, fontsize=14)
    map.savefig(map_path, format='svg', bbox_inches="tight")
    plt.close()

def process(base_path:str, data_path:str, output_path:str, heatmap_title: str, motifs_path: str):
    data_path=path.join(base_path,data_path)
    data = pd.read_csv(data_path)

    output_path=path.join(base_path,output_path)
    
    title=" ".join(heatmap_title.split(','))
    colors=[]
    if motifs_path:
        file_motifs=open(motifs_path)
        motifs=[]
        for line in file_motifs:
            items=line.split()
            motifs.append(items[0])
            colors.append(items[1])

        columns = ['sample_name'] + motifs
        data=data[columns]
    colors_list=colors if motifs_path else ''
    generate_heatmap(output_path, data, title, colors_list)


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('base_path', type=str, help='A path in which there is a subfolder for each bc')
    parser.add_argument('data_path',type=str,help='A path to csv file that conation motids and samples data')
    parser.add_argument('output_path', type=str, help='A path for create the heatmap')
    parser.add_argument('title_heatmap',type=str, help='String of the title heatmap, Separate each word with a comma: EX8,croos,best,motifs')
    parser.add_argument('--motifs', type=str,help='A path for motifs file and color')
    args = parser.parse_args()

    process(args.base_path,args.data_path,args.output_path,args.title_heatmap,args.motifs)

