import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
import sys
import os
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from tools.generate_heat_map import *
from auxiliaries.pipeline_auxiliaries import *

color_names=['red','green','blue','yellow','grey','pink','black','orange','purple','gold','brown']
def unit_bc_data(sample2bc_path, path_output, path_model_fitting, is_hits, num_of_motifs_evey_bc,selected_motifs,verbose):    
    #get the sample names form the file sample2bc
    sample2bc=load_table_to_dict(sample2bc_path,'More then one sample with the same name')
    sample_name=sorted(sample2bc.keys())
    bc=list(set(sample2bc.values()))    
    print(is_hits)
    dic_motif={}
    if selected_motifs:
        file_motif=pd.read_csv(selected_motifs)
        #transpose the data frame to a dictionary
        dic_motif=file_motif.set_index('name_bio').T.to_dict('list')
        color_list=False
    else:
        #color = list(np.random.choice(range(0, 10), size=len(bc)))
        #color=[(x*0.1,x*0.1,x*0.1) for x in color]
        color_list=[]
        for num,name_bc in enumerate(bc):
            path_motifs=path_model_fitting+'/'+name_bc+'/'+name_bc+'_hits_model/single_feature_accuracy.txt' if is_hits else path_model_fitting+'/'+name_bc+'/'+name_bc+'_values_model/single_feature_accuracy.txt'
            motifs=[]
            with open(path_motifs) as myfile:
                for line in myfile:
                    key,value= line.split('\t')
                    value_split=value.split('\n')[0]
                    if key=='Feature':
                        continue
                    if value_split!='1.0':
                        break
                    motifs.append(key)
            size=0
            if num_of_motifs_evey_bc and len(motifs)>num_of_motifs_evey_bc:
                size=len(motifs[:num_of_motifs_evey_bc])
                if size!=0:
                    dic_motif[name_bc]=motifs[:num_of_motifs_evey_bc]
            else:    
                size=len(motifs)
                if size!=0:
                    dic_motif[name_bc]=motifs
            for i in range(size):
                color_list.append(color_names[num])
    
    #get the data from every table of other biological condition
    df=pd.DataFrame({'sample_name':sample_name})
    for key in dic_motif.keys():
        path_file=path_model_fitting+'/'+key+'/'+key+'_hits_model/perfect_feature_names.csv' if is_hits else path_model_fitting+'/'+key+'/'+key+'_values_model/perfect_feature_names.csv'
        if os.path.exists(path_file):
            table=pd.read_csv(path_file)
            list_motifs=[x for x in dic_motif[key] if str(x) != 'nan']
            table_new= table[list_motifs]
            df=pd.merge(df,table_new,right_index=True,left_index=True)

    title_map='Heatmap of best motifs from every biological condition'
    generate_heatmap(path_output,df,title_map,color_list)  
     

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('samplename2bc', type=str, help='A file of samplename2biologicalcondition.txt')
    parser.add_argument('output_heatmap', type=str, help='A output folder to put the heatmap result')
    parser.add_argument('model_fitting_results', type=str, help='folder of result of random forest run')
    parser.add_argument('--isHit', action='store_true', help='If to connect all the hits data or values data')
    parser.add_argument('--num_of_motifs_evey_bc', type=int, help='How many motifs to choose from every BC.')
    parser.add_argument('--file_of_select_motif', type=str, help='csv file with specific selected motif to connect together')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
    unit_bc_data(args.samplename2bc, args.output_heatmap, args.model_fitting_results,True if args.isHit else False, args.num_of_motifs_evey_bc, args.file_of_select_motif, True if args.verbose else False)