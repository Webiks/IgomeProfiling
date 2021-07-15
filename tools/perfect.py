import pandas as pd

def simple_separation(file_path: str, filterIn: str = None, filterOut: str = None, bc_override: str = None):
    data = pd.read_csv(file_path, engine='python')
    if filterIn:
        data = data[data.sample_name.str.contains(filterIn)]
    if filterOut:
        data = data[~data.sample_name.str.contains(filterOut)]
    if bc_override:
        data.loc[data.sample_name.str.contains(bc_override), 'label'] = 'bc'
    data = data.transpose()
    # print(data['sample_name'][0])
    j = 0
    labels = []
    bc_label = ''
    perfect_motifs = []
    artifact_motifs = []
    mixed_motifs = []
    for tup in data.itertuples():
        j += 1
        if j == 1: # Samples
            print('labels', tup[1:])
            continue
        if j == 2: # groups 
            labels = list(tup[1:])
            groups = pd.unique(labels)
            print('groups', groups)
            bc_label = [group for group in groups if group != other_label][0]
            continue
        motif = tup[0]
        values = { key: [] for key in groups }
        for index, value in enumerate(tup[1:]):
            values[labels[index]].append(value)
        
        bc_min = min(values[bc_label])
        bc_max = max(values[bc_label])
        other_min = min(values[other_label])
        other_max = max(values[other_label])

        if bc_min > other_max:
            perfect_motifs.append(motif)
        elif bc_max < other_min: # should be <= ?
            artifact_motifs.append(motif) 
        else:
            mixed_motifs.append(motif)
    
    return { 'perfect': perfect_motifs, 'artifact': artifact_motifs, 'mixed': mixed_motifs }


def separation_containment(values, hits):
    contained = { key: [] for key in values }
    for key in values:
        for motif in values[key]:
            if motif in hits[key]:
                contained[key].append(motif)
    return contained


values_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting/ferret_texas/ferret_texas_values.csv'
hits_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_fitting/ferret_texas/ferret_texas_hits.csv'
other_label = 'other'

filterIn = 'Sanofi' # None, 'mAb', 'IgG'
filterOut = None # 'Sanofi'
bc_override = None # 'Texas_2012' # 'cov001', 'cov002'
values_separation = simple_separation(values_path, filterIn, filterOut, bc_override)
hits_separation = simple_separation(hits_path, filterIn, filterOut, bc_override)
contained = separation_containment(values_separation, hits_separation)

print('perfects', values_separation['perfect'])
print('artifacts', values_separation['artifact'])
if len(values_separation['perfect']) == len(contained['perfect']):
    print('All perfects are contained in hits')
else:
    print('Contained in hits (perfect)', contained['perfect'])
    print('Not contained in hits (perfect)', [x for x in values_separation['perfect'] if x not in contained['perfect']])
