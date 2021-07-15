'''
Merge meme files from different samples/biological conditions
'''
from os import path, listdir, popen


def is_valid_dir(dir, prefixes, ignore_suffixes):
    is_valid_prefix = False
    for prefix in prefixes:
        if dir.startswith(prefix):
            is_valid_prefix = True
            break
    if not is_valid_prefix: return False

    is_valid_suffix = True
    for suffix in ignore_suffixes:
        if dir.endswith(suffix):
            is_valid_suffix = False
            break

    return is_valid_suffix


def merge(base_path, output_path, prefix, ignore_suffix):
    print('Merging...')
    is_first = True
    files = [path.join(base_path, dir, 'meme.txt') for dir in listdir(base_path) if \
        is_valid_dir(dir, prefix, ignore_suffix)]
    with open(output_path, 'w') as w:
        for file in files:
            with open(file, 'r') as r:
                if is_first:
                    data = r.readlines()  # First time include headers
                    is_first = False
                else:
                    data = r.readlines()[6:]  # Ignore meme header
            w.writelines(data)


def unite(meme_path, cluster_output_path):
    print('Uniting...')
    unite_pssm_script_path = './UnitePSSMs/UnitePSSMs'
    aln_cutoff = 20
    pcc_cutoff = 0.6

    cmd = f'{unite_pssm_script_path} -pssm {meme_path} -out {cluster_output_path} ' \
          f'-aln_cutoff {aln_cutoff} -pcc_cutoff {pcc_cutoff}'
    stream = popen(cmd)
    output = stream.read()
    print(output)

if __name__ == '__main__':
    base_path = '/mnt/d/workspace/data/webiks/igome/results/hiv/exp_7_8_top2_mi'
    prefix = ['Top_EC']
    ignore_suffix = ['_1', '_2', '_3']
    meme_output_path = path.join(base_path, 'merged_meme.txt')
    cluster_output_path = path.join(base_path, 'merged_clusters.csv')
    merge(base_path, meme_output_path, prefix, ignore_suffix)
    unite(meme_output_path, cluster_output_path)
