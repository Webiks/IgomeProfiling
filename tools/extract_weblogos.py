from os import path
from motif_inference.generate_weblogo import generate_weblogo as create_weblogo

memes_path = '/mnt/d/workspace/data/webiks/igome/results/exp10_covid/cov_2230_mAb_memes.txt'
cleaned_path = '/mnt/d/workspace/data/webiks/igome/results/exp10_covid/2230_cleaned/'
output_path = '/mnt/d/workspace/data/webiks/igome/results/exp10_covid/weblogos'
motifs = ['GGSLRPGSS', 'GGSLRPTASS', 'CGSLKGTAQAPC', 'CHSMADITKGGC', 'CGSLRLSSLAF', 'CGSLKLTPSNSLY', 'CGRLDPFLNC', 'CGSLRACT', 'PYRKNPVDGLTP']

def extract_mapping(file_path: str):
    mapping = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('MOTIF '):
            parts = line.split(' ')[1].split('_')
            mapping[parts[0]] = '_'.join(parts[1:]).rstrip()
    return mapping

def generate_weblogo(motif: str):
    # Extract meme name
    faa_file = mapping[motif]
    meme_path = path.join(cleaned_path, faa_file)
    weblogo_output_path = path.join(output_path, f'MOTIF_{motif}_{faa_file}')
    print(f'Generating weblogo for {motif}')
    create_weblogo(meme_path, weblogo_output_path, 'NO TITLE')

mapping = extract_mapping(memes_path)
for motif in motifs:
    generate_weblogo(motif)
