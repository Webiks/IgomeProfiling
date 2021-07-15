import time
import sys
import os
import shutil
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')

if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import submit_pipeline_step, wait_for_results
from auxiliaries.pipeline_auxiliaries import log_scale


def parse_data(file_path):
    try:
        data = pd.read_csv(file_path, engine='python')
    except Exception as e:
        exit(f'Cannot analyse data! {e}')
    
    # create masks for separating train and test samples (fifth sample of each mAb)    
    test_rows_mask = data['sample_name'].str.contains('test')
    train_rows_mask = ~test_rows_mask

    # separating train and test samples (fifth sample of each mAb)
    # a matrix of the actual feature values
    y_train = data['label'][train_rows_mask]
    y_test = data['label'][test_rows_mask]

    sample_names_train = data['sample_name'][train_rows_mask]
    sample_names_test = data['sample_name'][test_rows_mask]

    # don't need the labels and sample names anymore
    data.drop(['sample_name', 'label'], axis=1, inplace=True)
    # a matrix of the actual feature values
    X_train = data[train_rows_mask].values
    X_test = data[test_rows_mask].values

    feature_names = np.array(data.columns)

    return X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test


def get_hyperparameters_grid(seed):
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start=100, stop=2000, num=20)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num=11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node
    min_samples_split = [2, 5, 10]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4, 8]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                   'max_features': max_features,
                   'max_depth': max_depth,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}
                   
    if os.path.exists('/Users/Oren'):
        # use all cores when running locally.
        # does not apply on the cluster (needed to be set as well in the .pbs file)
        random_grid['n_jobs'] = [-1]

    # Use the random grid to search for best hyperparameters
    return random_grid


def sample_configurations(hyperparameters_grid, num_of_configurations_to_sample, seed):
    configurations = []
    for i in range(num_of_configurations_to_sample):
        configuration = {}
        for key in hyperparameters_grid:
            np.random.seed(seed+i)
            configuration[key] = np.random.choice(hyperparameters_grid[key], size=1)[0]
        configurations.append(configuration)
    return configurations



def generate_heat_map(df, number_of_features, rank_method, number_of_samples, output_path):
    
    train_data = log_scale(df, rank_method)
    cm = sns.clustermap(train_data, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top {number_of_features} discriminatory motifs")
    cm.savefig(f"{output_path}.svg", format='svg', bbox_inches="tight")
    plt.close()



def save_model_features(X, feature_indexes, feature_names, sample_names, output_path):
    df = pd.DataFrame()
    for i in range(len(feature_names)):
        df[feature_names[i]] = X[:, feature_indexes[i]]  # train_data.iloc[:, :number_of_features]
    df.set_index(sample_names, inplace=True)
    df.to_csv(f"{output_path}.csv")
    return df



def write_results_feature_selection_summary(feature_selection_summary_path, path_dir):
    #feature_selection_summary_f = open(feature_selection_summary_path, 'a')
    models = sorted([x[0] for x in os.walk(path_dir)])
    del models[0]
    with open(feature_selection_summary_path, 'a') as feature_selection_summary_f:
        for path_number_model in models:        
            path_file = f'{path_number_model}/feature_selection.txt'
            if os.path.exists(path_file):
                with open(path_file) as infile:
                    line = infile.readline()
                    feature_selection_summary_f.write(line)
                os.remove(path_file)
    

def train_models(csv_file_path, done_path, logs_dir,error_path, num_of_configurations_to_sample, number_parallel_random_forest, min_value_error,
                 rank_method, cv_num_of_splits, seed, random_forest_seed, queue_name, verbose, argv):
    logging.info('Preparing output path...')
    csv_folder, csv_file_name = os.path.split(csv_file_path)
    csv_file_prefix = os.path.splitext(csv_file_name)[0]  # without extension
    output_path = os.path.join(csv_folder, f'{csv_file_prefix}_model')
    os.makedirs(output_path, exist_ok=True)

    best_model_path = os.path.join(output_path, f'best_model')

    feature_selection_summary_path = f'{output_path}/feature_selection_summary.txt'

    logging.info('Parsing data...')

    X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test = parse_data(csv_file_path)

    # single feature analysis
    logging.info('Applying single feature analysis...')
    perfect_feature_names, perfect_feature_indexes = measure_each_feature_accuracy(X_train, y_train, feature_names, output_path, seed, cv_num_of_splits)
    if perfect_feature_names:
        df = save_model_features(X_train, perfect_feature_indexes, perfect_feature_names, sample_names_train, f'{output_path}/perfect_feature_names')
        generate_heat_map(df, df.shape[1], rank_method, df.shape[0], f'{output_path}/perfect_feature_names')
    else:
        # touch a file so we can see that there were no perfect features
        with open(f'{output_path}/perfect_feature_names', 'w') as f:
            pass

    # feature selection analysis
    logging.info('\nApplying feature selection analysis...')
    if cv_num_of_splits < 2:
        logging.info('Number of CV folds is less than 2. '
                     'Updating number of splits to number of samples and applying Leave One Out approach!')
        cv_num_of_splits = len(y_train)
    logging.info(f'Number of CV folds is {cv_num_of_splits}.')

    logger.info('\n'+'#'*100 + f'\nTrue labels:\n{y_train.tolist()}\n' + '#'*100 + '\n')

    logging.info('Sampling hyperparameters...')
    hyperparameters_grid = get_hyperparameters_grid(seed)

    feature_selection_summary_f = open(feature_selection_summary_path, 'w')
    feature_selection_summary_f.write(f'model_number\tnum_of_features\tfinal_error_rate\n')
    feature_selection_summary_f.close()

    sampled_configurations = sample_configurations(hyperparameters_grid, num_of_configurations_to_sample, seed)
    script_name = 'model_fitting/train_random_forest.py'
    num_of_expected_results = 0
    all_cmds_params = []
    for i, configuration in enumerate(sampled_configurations):
        model_number = str(i).zfill(len(str(num_of_configurations_to_sample)))
        output_path_i = os.path.join(output_path, model_number)
        logging.info(f'Creating output path #{i}...')
        os.makedirs(output_path_i, exist_ok=True)
        file_save_configuration = save_configuration_to_txt_file(configuration, output_path_i)
        logging.info(f'Configuration #{i} hyper-parameters are:\n{configuration}')
        done_file_path=os.path.join(logs_dir, f'{model_number}_{csv_file_prefix}_done_train_random_forest.txt')
        cmds=[file_save_configuration, csv_file_path, rank_method, output_path_i,
                model_number, done_file_path, '--random_forest_seed', random_forest_seed,
                '--cv_num_of_splits', cv_num_of_splits]
        if not os.path.exists(done_file_path):
            all_cmds_params.append(cmds)
        else: 
            logger.debug(f'Skipping random forest train as {done_file_path} found')
            num_of_expected_results +=1

    executable = 'python'  
    if len(all_cmds_params) > 0:
        stop = False
        for count, cmds_params in enumerate(all_cmds_params):
            model_number = cmds_params[4]
            cmd = submit_pipeline_step(script_name, [cmds_params],
                                logs_dir, f'{csv_file_prefix}_{model_number}_train_random_forest',
                                queue_name, verbose, executable=executable)
            num_of_expected_results += 1  # a single job for each train random forest
            if (count + 1) % number_parallel_random_forest == 0 or (count + 1) == num_of_configurations_to_sample:
                wait_for_results(script_name, logs_dir, num_of_expected_results, example_cmd=cmd,
                        error_file_path = error_path, suffix = f'_{csv_file_prefix}_done_train_random_forest.txt') 
                #write the results to feature_selection_summary file
                write_results_feature_selection_summary(feature_selection_summary_path, output_path)
                #check if we found the best model and can stop run
                models_stats = pd.read_csv(feature_selection_summary_path, sep='\t', dtype={'model_number': str, 'num_of_features':int, 'final_error_rate': float })
                for feature, error in zip(models_stats['num_of_features'], models_stats['final_error_rate']):
                    if feature == 1 and error <= min_value_error:
                        stop = True
                        break
                if stop:
                    break
    else:
         logger.info(f'Skipping random forest train, all found')       

    feature_selection_summary_f.close()

    # find who was the best performing model
    models_stats = pd.read_csv(feature_selection_summary_path, sep='\t', dtype={'model_number': str, 'num_of_features':int, 'final_error_rate': float })
    lowest_error_models = models_stats[models_stats['final_error_rate'] ==
                                       min(models_stats['final_error_rate'])]

    # in case there is a tie between best models, we choose the model with the lowest number of features
    if len(lowest_error_models) > 1:
        lowest_error_models = lowest_error_models[lowest_error_models['num_of_features'] ==
                                                  min(lowest_error_models['num_of_features'])]

    best_model = lowest_error_models['model_number'].iloc[0]

    # keep only best model (and remove the rest)
    shutil.copytree(f'{output_path}/{best_model}', best_model_path)
    for folder in os.listdir(output_path):
        path = f'{output_path}/{folder}'
        if folder == 'best_model' or not os.path.isdir(path):
            continue
        shutil.rmtree(path)

    with open(done_path, 'w') as f:
        f.write(' '.join(argv) + '\n')


def measure_each_feature_accuracy(X_train, y_train, feature_names, output_path, seed, cv_num_of_splits):

    feature_to_avg_accuracy = {}    
    rf = RandomForestClassifier(random_state=np.random.seed(seed))

    for i, feature in enumerate(feature_names):
        logger.info(f'Checking feature {feature} number {i}')
        cv_score = cross_val_score(rf, X_train[:, i].reshape(-1, 1), y_train, cv=StratifiedKFold(cv_num_of_splits, shuffle=True)).mean()
        if cv_score == 1:
            logger.info('-' * 10 + f'{feature} has 100% accuracy!' + '-' * 10)
        feature_to_avg_accuracy[feature] = cv_score

    with open(f'{output_path}/single_feature_accuracy.txt', 'w') as f:
        f.write('Feature\tAccuracy_on_cv\n')
        for feature in sorted(feature_to_avg_accuracy, key=feature_to_avg_accuracy.get, reverse=True):
            f.write(f'{feature}\t{feature_to_avg_accuracy[feature]}\n')

    perfect_feature_names = []
    perfect_feature_indexes = []
    with open(f'{output_path}/features_with_perfect_accuracy.txt', 'w') as f:
        for i, feature in enumerate(feature_to_avg_accuracy):
            if feature_to_avg_accuracy[feature] == 1:
                perfect_feature_names.append(feature)
                perfect_feature_indexes.append(i)
                f.write(f'{feature}\n')

    return perfect_feature_names, perfect_feature_indexes


def save_configuration_to_txt_file(sampled_configuration, output_path_i):
    file_save_configuration=f'{output_path_i}/hyperparameters_configuration.txt'
    with open(f'{output_path_i}/hyperparameters_configuration.txt', 'w') as f:
        for key in sampled_configuration:
            f.write(f'{key}={sampled_configuration[key]}\n')
    joblib.dump(sampled_configuration, f'{output_path_i}/hyperparameters_configuration.pkl')
    return file_save_configuration

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model ')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('logs_dir', help='A path for the log dir')
    parser.add_argument('error_path', help='Path for error file')
    
    parser.add_argument('--num_of_configurations_to_sample', default=100, type=int, help='How many random configurations of hyperparameters should be sampled?')
    parser.add_argument('--number_parallel_random_forest', default=20, type=int, help='How many random forest configurations to run in parallel')
    parser.add_argument('--min_value_error_random_forest', default=0, type=float, help='A random forest model error value for convergence allowing to stop early')
    parser.add_argument('--cv_num_of_splits', type=int ,default=2, help='How folds should be in the cross validation process? (use 0 for leave one out)')
    parser.add_argument('--seed', type=int, default=42, help='Seed number for reconstructing experiments')    
    parser.add_argument('--random_forest_seed', default=123 , type=int, help='Random seed value for generating random forest configurations')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles', 'hits'], default='hits', help='Motifs ranking method')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    train_models(args.data_path, args.done_file_path, args.logs_dir, args.error_path, args.num_of_configurations_to_sample, args.number_parallel_random_forest,
                 args.min_value_error_random_forest, args.rank_method, args.cv_num_of_splits, args.seed, args.random_forest_seed, args.queue ,args.verbose, argv=sys.argv)
