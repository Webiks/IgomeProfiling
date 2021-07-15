import re

cluster_pattern = re.compile('(.+?)_(.+)_clusterRank')

def test(path: str):
    with open(path, 'r') as f:
        lines = f.readlines()
    i = 0
    for line in lines:
        i += 1
        clusters = line.split(',')
        biological_condition_groupping = {}
        for cluster in clusters:
            cluster_match = cluster_pattern.match(cluster)
            consensus = cluster_match.group(1)
            biological_condition = cluster_match.group(2)
            try:
                biological_condition_groupping[biological_condition].append(consensus)
            except:
                biological_condition_groupping[biological_condition] = [consensus]
        if len(biological_condition_groupping.keys()) > 1:
            print(f'Cluster {i}, {biological_condition_groupping}')


if __name__ == '__main__':
    file_path = '/mnt/d/workspace/data/webiks/igome/results/exp4+9_all_half_gaps/model_inference/cross/texas_combined.csv'
    test(file_path)
