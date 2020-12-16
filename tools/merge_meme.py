'''
Merge meme files from different samples/biological conditions
'''
import sys
from os import path, listdir, popen
from auxiliaries.pipeline_auxiliaries import *

def is_valid_dir(dir,prefixes):
    is_valid_prefix = False
    for prefix in prefixes:
        if dir.startswith(prefix):
            is_valid_prefix = True
            break
    return is_valid_prefix        

def merge(base_path,s2b_path,output_path,prefixes):
    print('Merging...')
    s2b=load_table_to_dict(s2b_path,'Barcode {} belongs to more than one sample!!')
    list_prefixes=prefixes.split(',')
    bc=list(set(s2b.values()))
    is_first = True
    if list_prefixes[0]=='all_bc':
        files = [path.join(base_path, name_bc, 'meme.txt') for name_bc in bc]        
    else:
        files = [path.join(base_path, name_bc, 'meme.txt') for name_bc in bc if \
            is_valid_dir(name_bc, list_prefixes)]
    with open(output_path, 'w') as w:
        for file in files:
            with open(file, 'r') as r:
                if is_first:
                    data = r.readlines()  # First time include headers
                    is_first = False
                else:
                    data = r.readlines()[6:]  # Ignore meme header
            w.writelines(data)


def unite(meme_path, cluster_output_path,aln_cutoff,pcc_cutoff,unite_pssm_script_path):
    print('Uniting...')
    cmd = f'{unite_pssm_script_path} -pssm {meme_path} -out {cluster_output_path} ' \
          f'-aln_cutoff {aln_cutoff} -pcc_cutoff {pcc_cutoff}'
    stream = popen(cmd)
    output = stream.read()
    print(output)

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('base_path',type=str ,help='A path for results of motif inference phase.')
    parser.add_argument('samples2bc',type=str,help='A path for txt file that contain the samples name and bc names.')
    parser.add_argument('--prefixes',type=str, default='all_bc',help='A string of prefix of bc to merge. should write with comma separatar.For exmple: Top_,EC_')
    parser.add_argument('unite_pssm_script_path',type=str,help='A path for unite_pssm ./UnitePSSMs/UnitePSSMs')
    parser.add_argument('--aln_cutoff',type=int,default=20)
    parser.add_argument('--pcc_cutoff',type=float,default=0.6)
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    meme_output_path = path.join(args.base_path, 'merged_meme.txt')
    cluster_output_path = path.join(args.base_path, 'merged_clusters.csv')
    merge(args.base_path,args.samples2bc ,meme_output_path,args.prefixes)
    unite(meme_output_path, cluster_output_path,args.aln_cutoff,args.pcc_cutoff,args.unite_pssm_script_path)
