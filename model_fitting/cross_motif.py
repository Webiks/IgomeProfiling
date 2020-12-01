import pandas as pd
import os
import sys

def cross_motif_results(merge_cluster_path,results_dir_path,done_file_path):
    #create a file for output
    f_out=open(done_file_path, 'w')
    merge_cluster=open(merge_cluster_path,'r')
    result={}
    #for every EC_HIV open the file distinctive_motifs_top_10.csv
    for dirname in os.listdir(results_dir_path):
        ### if only between 10 motifs
        #result={}
        df=pd.read_csv(f'{results_dir_path}/{dirname}/{dirname}_values_distinctive_motifs/distinctive_motifs_top_10.csv')
        motifs_name=df['motif']
        # find witch motif in this file are connected for the same cluster.
        for num,line in enumerate(merge_cluster):
            line_var=line.split(',')
            for var in line_var:
                motif=var.split('_')[0]
                for name in motifs_name:
                    if name==motif:
                        new_name_motif=motif +' '+dirname
                        if num in result.keys() and new_name_motif not in result[num]:
                            result[num].append(new_name_motif)
                        else:
                            result[num]=[new_name_motif]
        merge_cluster.seek(0, 0)
    #f_out.write(f'Name folder: {dirname')
    #print result to the file
    print(result)
    for key in result.keys():
        if len(result[key])>1:
            f_out.write(f'number line: {key} motifs: {result[key]}\n')
    f_out.close()
    merge_cluster.close()




if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('merge_cluster_path', help='A path to a csv file with cluster for meme that merge all the groups.')
    parser.add_argument('results_dir_path', help='A path in which each csv file contains hit/value distincte motif. ')
    parser.add_argument('done_file_path', help='A path to a file contain the results.')
    args = parser.parse_args()

    cross_motif_results(args.merge_cluster_path, args.results_dir_path,args.done_file_path)

