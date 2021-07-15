import json 
import sys
import seaborn as sns
import joblib
import os
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import plot_roc_curve
import pandas as pd
import numpy as np
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)
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


def generate_roc_curve(X, y, classifier, number_of_features, output_path):
    ax = plt.gca()
    plot_roc_curve(classifier, X, y, **{'marker': 'o'}, ax=ax)
    plt.savefig(f"{output_path}/roc_curve_{number_of_features}.png", format='png', bbox_inches="tight")
    plt.close()


def generate_heat_map(df, number_of_features, rank_method, number_of_samples, output_path):
    train_data = log_scale(df, rank_method)
    cm = sns.clustermap(train_data, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top {number_of_features} discriminatory motifs")
    cm.savefig(f"{output_path}.svg", format='svg', bbox_inches="tight")
    plt.close()

def plot_error_rate(errors, features, cv_num_of_splits, output_path_dir):
    raw_data_path = f"{output_path_dir}/error_rate.txt"
    with open(raw_data_path, 'w') as f:
        f.write('Features\tErrors\n')
        f.write('\n'.join(f'{feature}\t{error}' for feature, error in zip(features[::-1], errors[::-1])))
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel('# of Features')
    plt.ylabel(f'{cv_num_of_splits}Fold CV Avg. Error Rate')
    plot_path = raw_data_path.replace('.txt', '.png')
    plt.savefig(plot_path)
    plt.close()

def train(rf, X, y, feature_names, sample_names, rank_method, output_path, cv_num_of_splits):

    original_feature_names = feature_names[:]
    original_X = X[:]
    logger.debug('\n'+'#'*100 + f'\nTrue labels:\n{y.tolist()}\n' + '#'*100 + '\n')

    # Fit the best model configuration on the WHOLE dataset
    logging.info('Training...')
    model = rf.fit(X, y)
    importance = model.feature_importances_

    # the permutation needed to get the feature importances in a decreasing order
    decreasing_feature_importance = np.argsort(importance)[::-1]
    assert (sorted(importance, reverse=True) == importance[decreasing_feature_importance]).all()

    # the indexes that will be used in the next fitting process
    # At first, we start with all features. Next we will remove less important once.
    features_indexes_to_keep = range(len(feature_names))

    # write feature importance to storage
    with open(f'{output_path}/feature_importance.txt', 'w') as f:
        for i in range(len(importance)):
            f.write(f'{feature_names[i]}\t{importance[i]}\n')

    # write sorted feature importance to storage
    with open(f'{output_path}/sorted_feature_importance.txt', 'w') as f:
        for i in range(len(importance)):
            f.write(f'{feature_names[decreasing_feature_importance[i]]}\t'
                    f'{importance[decreasing_feature_importance[i]]}\n')

    number_of_samples, number_of_features = X.shape
    cv_avg_error_rate = previous_cv_avg_error_rate = 1
    cv_avg_error_rates = []
    number_of_features_per_model = []
    while cv_avg_error_rate <= previous_cv_avg_error_rate and number_of_features >= 1:

        # save previous cv_avg_error_rate to make sure the performances do not deteriorate
        previous_cv_avg_error_rate = cv_avg_error_rate

        # compute current model accuracy for each fold of the cross validation
        cv_score = cross_val_score(model, X, y, cv=StratifiedKFold(cv_num_of_splits))

        # current model cv_avg_error_rate rate
        cv_avg_error_rate = 1 - cv_score.mean()
        number_of_features_per_model.append(number_of_features)
        cv_avg_error_rates.append(cv_avg_error_rate)

        logger.info(f'Number of features is {number_of_features} with avg. error rate of {cv_avg_error_rate}')
        if cv_avg_error_rate > previous_cv_avg_error_rate:
            # Stop procedure and discard current stats
            break

        # save current model (unsorted) features to a csv file
        df = save_model_features(original_X, features_indexes_to_keep, feature_names, sample_names, f'{output_path}/Top_{number_of_features}_features')

        generate_heat_map(df, number_of_features, rank_method, number_of_samples, f'{output_path}/{number_of_features}')

        generate_roc_curve(X, y, rf, number_of_features, output_path)

        # save the model itself (serialized) for future use
        joblib.dump(model, os.path.join(output_path, f'Top_{number_of_features}_features_model.pkl'))

        # Sanity check for debugging: predicting the test data
        model_score = model.score(X, y)
        if 1-model_score > cv_avg_error_rate:
            predictions = model.predict(X).tolist()
            logging.error('1-model_score > cv_avg_error_rate !!!')
            logger.info(f'Full model error rate is {1-model_score}')
            logger.info(f'Current model\'s predictions\n{predictions}')
            logger.info(f'number_of_features {number_of_features}')
            logger.info(f'output_path {output_path}')

        # update number of features
        number_of_features //= 2

        # extract only the (new) half most important features
        features_indexes_to_keep = sorted(decreasing_feature_importance[:number_of_features])
        feature_names = original_feature_names[features_indexes_to_keep]
        X = original_X[:, features_indexes_to_keep]

        if number_of_features > 0:
            # re-evaluate
            model = rf.fit(X, y)

    return cv_avg_error_rates, number_of_features_per_model
def save_model_features(X, feature_indexes, feature_names, sample_names, output_path):
    df = pd.DataFrame()
    for i in range(len(feature_names)):
        df[feature_names[i]] = X[:, feature_indexes[i]]  # train_data.iloc[:, :number_of_features]
    df.set_index(sample_names, inplace=True)
    df.to_csv(f"{output_path}.csv")
    return df

def configuration_from_txt_to_dictionary(configuration_path):
    configuration = {}
    with open(configuration_path, 'r') as f:
            for line in f:
                (key, val) = line.split('=')
                if key == 'max_features' or key == 'bootstrap':
                    configuration[key] = val.split('\n')[0]
                elif (key == 'max_depth' and val == 'None\n'):
                    configuration[key]=None
                else:
                    configuration[key] = int(val)
    return configuration                

def pre_train(configuration_path, csv_file_path, rank_method, output_path_i, model_number, done_file_path, random_forest_seed, cv_num_of_splits, argv):
    configuration = configuration_from_txt_to_dictionary(configuration_path)
    logger.info(f'Start run random forest for model number {model_number}:')
    feature_selection_f = open(f'{output_path_i}/feature_selection.txt', 'w')
    rf = RandomForestClassifier(**configuration, random_state=random_forest_seed)
    X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test = parse_data(csv_file_path)
    errors, features = train(rf, X_train, y_train, feature_names, sample_names_train, rank_method, output_path_i, cv_num_of_splits)
    plot_error_rate(errors, features, cv_num_of_splits, output_path_i)
    feature_selection_f.write(f'{model_number}\t{features[-1]}\t{errors[-1]}\n')
    feature_selection_f.close()
    
    with open(done_file_path, 'w') as f:
        f.write(' '.join(argv) + '\n')

if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('configuration_path', type=str, help='A dictionary of the configuration to this model ')
    parser.add_argument('csv_file_path', type=str, help='A csv file with data matrix to model')
    parser.add_argument('rank_method', type=str, help='if is hits data or pvalue data')
    parser.add_argument('output_path_i', type=str, help='A path for the results of this model')
    parser.add_argument('model_number', type=str, help='number of model')
    parser.add_argument('done_file_path', type=str, help='A path to a file that signals that the script finished running successfully')
    parser.add_argument('--random_forest_seed', default=123 , type=int, help='Random seed value for generating random forest configurations')
    parser.add_argument('--cv_num_of_splits', default=2, type=int, help='number of CV folds')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()
    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

   
    pre_train(args.configuration_path, args.csv_file_path, args.rank_method, args.output_path_i, args.model_number, args.done_file_path,
             args.random_forest_seed, args.cv_num_of_splits, argv=sys.argv)
