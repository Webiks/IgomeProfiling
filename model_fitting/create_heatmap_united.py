import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
import sys
import os


def getData(sample2bc,path_output,path_model_fitting,selected_motifs,isHits):    
    dic_motif={}
    if selected_motifs:
        file_motif=pd.read_csv(selected_motifs)
        #transpose the data frame to a dictionary
        dic_motif=file_motif.set_index('name_bio').T.to_dict('list')
    #get the sample names form the file sample2bc
    sample2bc_file=open(sample2bc,'r')
    sample_name=[]
    bc=[]
    for var in sample2bc_file:
        var_split=var.split()
        sample_name.append(var_split[0])
        bc.append(var_split[1])
    bc=list(set(bc))
    #get the hits/value from every table of other biological condition
    data=pd.DataFrame({'sample_name':sample_name})
    for name_bc in bc:
        print(name_bc)
        path_file=''
        if isHits:
            path_file=path_model_fitting+'/'+name_bc+'/'+name_bc+'_hits_model/perfect_feature_names.csv'
        else:
            path_file=path_model_fitting+'/'+name_bc+'/'+name_bc+'_values_model/perfect_feature_names.csv'
        if os.path.exists(path_file):
            table=pd.read_csv(path_file)
            if not dic_motif:
                #the dictionary is empty, take all the motifs
                data=pd.merge(data,table)
            else:
                values=dic_motif[key]
                clean_values=[x for x in values if str(x) != 'nan']
                table_new= table[clean_values]
                data=pd.merge(data,table_new)
        print(data)
    motif_name=data.columns
    create_heatmap(data,sample_name,motif_name,path_output,isHits)
    
def create_heatmap(df,sample_name,motif_name,path_output,isHits):
    print(df)
    if isHits:
        #normalization the data, the data is mix of many runs. 
        df= df.drop(columns=['sample_name'])
        train_data=np.log2(df+1)
        train_data.insert(loc=0, column='sample_name', value=sample_name)
    else:
        train_data=df
    if not os.path.exists(path_output):
        os.makedirs(path_output, exist_ok=True)
    train_data.to_csv(path_output+'/data.csv',index=False)
    df_new=pd.read_csv(path_output+'/data.csv',index_col=0)
    cm = sb.clustermap(df_new, cmap="Blues",col_cluster=False)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top motifs")
    cm.savefig(f"{path_output}.svg", format='svg', bbox_inches="tight")
    plt.close()
    
    
     

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('samplename2bc', type=str, help='A file of samplename2biologicalcondition.txt')
    parser.add_argument('output_heatmap', type=str, help='A output folder to put the heatmap result')
    parser.add_argument('model_fitting_results', type=str, help='folder of result of random forest run')
    parser.add_argument('--isHits', default=True ,action='store_true')
    parser.add_argument('--file_of_select_motif', type=str, help='csv file with specific selected motif to connect together')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')
    
    getData(args.samplename2bc,args.output_heatmap,args.model_fitting_results,args.file_of_select_motif,True if args.isHits else False)
