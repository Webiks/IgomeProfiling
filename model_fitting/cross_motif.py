import pandas as pd
import os
import sys

def cross_motif_results(merge_cluster_path,results_dir_path,done_file_path,file_path_for_heatmap):
    #create a file for output
    f_out=open(done_file_path, 'w')
    f_out_for_heatmap=open(file_path_for_heatmap,'w')
    merge_cluster=open(merge_cluster_path,'r')
    result={}
    merge_cluster=open(merge_cluster_path,'r')
    #for every bc open the file distinctive_motifs_top_10.csv
    for dirname in os.listdir(results_dir_path):
        merge_cluster=open(merge_cluster_path,'r')
        df=pd.read_csv(f'{results_dir_path}/{dirname}/{dirname}_values_distinctive_motifs/distinctive_motifs_top_10.csv')
        motifs_name=df['motif']
        # find witch motif in this file are connected for the same cluster.
        for num,line in enumerate(merge_cluster,1):
            line_var=line.split(',')
            for var in line_var:
                motif=var.split('_')[0]
                for name in motifs_name:
                    if name==motif:
                        new_name_motif=motif +' '+dirname
                        if num in result.keys():
                            if new_name_motif in result[num]:
                                continue
                            else:    
                                result[num].append(new_name_motif)
                        else:
                            result[num]=[new_name_motif]
        merge_cluster.seek(0)
    for key in result.keys():
        if len(result[key])>1:
            f_out.write(f'number line: {key} motifs: {result[key]}\n')
            var=(result[key])[0].split()
            #from every cluster same the firsy name motifs to create heatmap
            f_out_for_heatmap.write(f'{var[1]} {var[0]}')
    f_out.close()
    f_out_for_heatmap.close()
    merge_cluster.close()



if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('merge_cluster_path', help='A path to a csv file with cluster for meme that merge all the groups.')
    parser.add_argument('results_dir_path', help='A path in which each csv file contains hit/value distincte motif. ')
    parser.add_argument('done_file_path', help='A path to a file contain the results.')
    parser.add_argument('file_path_for_heatmap', help='A path to a file contain motifs for create heatmap.')
    args = parser.parse_args()

    cross_motif_results(args.merge_cluster_path, args.results_dir_path,args.done_file_path,file_path_for_heatmap)

    #python3 model_fitting/cross_motif.py output/analysis/merged_clusters.csv output/analysis/model_fitting output/analysis/result_file.txt output/analysis/motifs_name_file.txt