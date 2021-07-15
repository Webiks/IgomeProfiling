import os
from sys import argv
first_phase_output_dir = argv[1]
path_to_samplename2biologicalcondition = argv[2]
resuls_path = argv[3]
result=''
for sample_name in os.listdir(first_phase_output_dir):
    samplename2biologicalcondition = {}
    with open(path_to_samplename2biologicalcondition) as f:
        for line in f:
            if line.strip():
                tokens = line.strip().split()
                samplename2biologicalcondition[tokens[0]] = tokens[1]

    sample_path=os.path.join(first_phase_output_dir, sample_name)
    print('file name is: ' + str(sample_name))
    if not os.path.isdir(sample_path):
        # print('Skipping sample_name ' + str(sample_name) + ' (not a directory).')
        continue
    file = ''
    file_name_is_legal = False
    for file_name in os.listdir(sample_path):
        if 'UNIQUE.NORMALIZATION_RpM' in file_name:
            file_name_is_legal = True
            break
    if not file_name_is_legal:
        print('Warning!!!!!\nsample_name ' + str(sample_name) + ' does not contain UNIQUE.NORMALIZATION_RpM file.')
    file_path = os.path.join(sample_path, file_name)
    print('file_path is: ' + str(file_path))
    #if len(sample_name)>8 and sample_name[8] == '_': #e.g., TTACAGCG_Israel_HCV_RNAplus_02_a_Random
    sample_name_without_barcode = sample_name[9:]
    if sample_name_without_barcode in samplename2biologicalcondition:
        biological_condition = samplename2biologicalcondition[sample_name_without_barcode]
        result += file_path + ',' + biological_condition + ',' + sample_name + '\n'
    else:
        print(f'{sample_name_without_barcode} is not in samplename2biologicalcondition. Skipping...')
with open(resuls_path, 'w') as f:
    f.write(result)

'''
import shutil, os
def copy_files(barcode2samplename, first_phase_output, analysis_output):
    with open(barcode2samplename) as f:
        text=[line.split()[0] for line in f.readlines() if line.strip()]
    # print(text)
    for dir in os.listdir(first_phase_output):
        sample_barcode=dir.split('_')[0]
        print(sample_barcode)
        if sample_barcode in text:
            shutil.copytree(os.path.join(first_phase_output, dir), os.path.join(analysis_output, dir))

analysis = '7_mAb_c108_vs_Patient_14'; copy_files('/groups/pupko/orenavr2/gershoni/Experiments/Exp_DP_4/analyses/'+analysis+'/metadata/barcode2samplename.txt', '/groups/pupko/orenavr2/gershoni/Experiments/Exp_DP_4/first_phase_output', '/groups/pupko/orenavr2/gershoni/Experiments/Exp_DP_4/analyses/'+analysis+'/first_phase_output')
'''

'''
import pandas as pd
import os
from functools import reduce
results_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_4/analyses/4_Ferrets_Vac_vs_Vac/results'
results_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/Avia_TB/results'
results_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_5/analyses/Analysis1/results'
dfs = [pd.read_csv(results_dir+'/'+csv) for csv in os.listdir(results_dir) if csv.endswith('csv')]
final_df = reduce(lambda left,right: pd.merge(left,right,on='Sample_Name'), dfs) #merge all dfs
final_df = final_df.reindex(list([a for a in final_df.columns if a != 'Sample_Name'])+['Sample_Name'], axis=1) #move sample_name to last column
final_df.to_csv(results_dir+'/../merged_df.csv', index=False)
'''

'''
top_k_features=[f'SICK.{cluster}' for cluster in ['Cluster_182.75_MEMBERS_4091.8321895832_Peptides_Rank_190', 'Cluster_321.22_MEMBERS_2951.71815731683_Peptides_Rank_220', 'Cluster_228.54_MEMBERS_3460.99020204534_Peptides_Rank_209', 'Cluster_247.82_MEMBERS_6917.04360534905_Peptides_Rank_124', 'Cluster_134.285_MEMBERS_16737.0976732608_Peptides_Rank_68', 'Cluster_141.178_MEMBERS_22295.3954140836_Peptides_Rank_54', 'Cluster_316.24_MEMBERS_1545.97103841776_Peptides_Rank_324', 'Cluster_26.537_MEMBERS_8159.78363436453_Peptides_Rank_111', 'Cluster_323.22_MEMBERS_2556.81443154997_Peptides_Rank_246', 'Cluster_207.111_MEMBERS_5094.90561974032_Peptides_Rank_158', 'Cluster_69.392_MEMBERS_10655.8195248553_Peptides_Rank_90']] + ['Sample_Name']
import pandas as pd
input_table = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/Avia_TB/merged_df.csv'
output_table = input_table[:input_table.find('.csv')] + '_top_k_features.csv'
df = pd.read_csv(input_table)
final_df = df[top_k_features]
final_df.to_csv(output_table, index=False)
'''

'''
csvs=['P5E2_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv', 'P5E3_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_14_Features.csv', 'P12E2_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_4_Features.csv', 'P12E3_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv', 'P4E2_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_4_Features.csv', 'P4E3_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_49_Features.csv', 'P18E2_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_3_Features.csv', 'P18E3_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_96_Features.csv', 'P15E3_against_Exp_DP_2.TESTING_PVals_significant_binary_labeled.Top_28_Features.csv']
import pandas as pd
import os
from functools import reduce
results_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/Exp_DP_2/analyses/TB_per_patient/random_forest'
dfs = [pd.read_csv(results_dir+'/'+csv, index_col=0).drop('label', 1) for csv in csvs]
final_df = reduce(lambda left,right: pd.merge(left,right,on='sample_name'), dfs) #merge all dfs
final_df = final_df.reindex(list([a for a in final_df.columns if a != 'sample_name'])+['sample_name'], axis=1) #move sample_name to last column
final_df.to_csv(results_dir+'/TB_motifs_merged_pvals_df.csv', index=False)

'''
