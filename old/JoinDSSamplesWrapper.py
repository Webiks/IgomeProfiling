import os, sys, subprocess

def fetch_join_ds_samples(correction = '50', gershoni_dir = '/groups/pupko/orenavr2/gershoni'):
    experiments_dir = os.path.join(gershoni_dir, 'Experiments')
    py_script = os.path.join(gershoni_dir, 'src/JoinDSSamples.py')
    domains_of_interest_dir = os.path.join(gershoni_dir, 'domains_of_interest')

    for exp in os.listdir(experiments_dir):
        # if 'DP' not in exp:
        #     continue
        experiment_dir = os.path.join(gershoni_dir, 'Experiments', exp)
        cmds_file = os.path.join(experiment_dir, 'JoinDSSample.cmds')
        cmds = ''
        for domains_of_interest_file in os.listdir(domains_of_interest_dir):
            domains_of_interest_path = os.path.join(domains_of_interest_dir, domains_of_interest_file)
            cmds += ' '.join(['python','-u', py_script, domains_of_interest_path, os.path.join(experiment_dir, 'first_phase_output'), correction]) + '!@#'
        #cmds += '\t' + 'JoinDSSample
        with open(cmds_file, 'w') as f:
            f.write(cmds)
        system_call = '/groups/pupko/orenavr2/src/cmds_fetcher.py {}'.format(cmds_file)
        print(system_call)
        subprocess.check_output(system_call, shell=True)

if __name__ == '__main__':
    print('Usage: python '+sys.argv[0]+' <?correction[50]> <?gershoni_dir[/groups/pupko/orenavr2/gershoni]')
    fetch_join_ds_samples(*sys.argv[1:])

