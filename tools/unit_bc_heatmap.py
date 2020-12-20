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

def unit_bc_data(sample2bc_path,path_output,path_model_fitting,is_hits,selected_motifs,verbose):    
    dic_motif={}
    if selected_motifs:
        file_motif=pd.read_csv(selected_motifs)
        #transpose the data frame to a dictionary
        dic_motif=file_motif.set_index('name_bio').T.to_dict('list')
    #get the sample names form the file sample2bc
    sample2bc=load_table_to_dict(sample2bc_path,'More then one sample with the same name')
    sample_name=sorted(sample2bc.keys())
    bc=list(set(sample2bc.values()))
    
    #get the data from every table of other biological condition
    df=pd.DataFrame({'sample_name':sample_name})
    for name_bc in bc:
        path_file=''
        if is_hits:
            path_file=path_model_fitting+'/'+name_bc+'/'+name_bc+'_hits_model/perfect_feature_names.csv'
        else:
            path_file=path_model_fitting+'/'+name_bc+'/'+name_bc+'_values_model/perfect_feature_names.csv'
        if os.path.exists(path_file):
            table=pd.read_csv(path_file)
            if not dic_motif:
                #the dictionary is empty, take all the motifs
                df=pd.merge(df,table)
            else:
                list_motifs=[x for x in dic_motif[name_bc] if str(x) != 'nan']
                table_new= table[list_motifs]
                df=pd.merge(df,table_new)
    title_map='Heatmap unit of all the biological condition'
    generate_heatmap(path_output,df,title_map)  
     

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('samplename2bc', type=str, help='A file of samplename2biologicalcondition.txt')
    parser.add_argument('output_heatmap', type=str, help='A output folder to put the heatmap result')
    parser.add_argument('model_fitting_results', type=str, help='folder of result of random forest run')
    parser.add_argument('isHit', type=bool, help='If to connect all the hits data or values data')
    parser.add_argument('--file_of_select_motif', type=str, help='csv file with specific selected motif to connect together')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
    
    unit_bc_data(args.samplename2bc,args.output_heatmap,args.model_fitting_results,args.isHit,args.file_of_select_motif,True if args.verbose else False)