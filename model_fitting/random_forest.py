import matplotlib
matplotlib.use('Agg')

import sys
import os
import shutil
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, cross_validate
from sklearn.metrics import plot_roc_curve

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def parse_data(file_path):
    # reading the CSV file if it's legal
    try:
        data = pd.read_csv(file_path, engine='python')
    except Exception as e:
        exit(f'Cannot analyse data! {e}')

    # separating train and test samples (fifth sample of each mAb)
    train_data = data[~data['sample_name'].str.contains('test')]
    test_data = data[data['sample_name'].str.contains('test')]

    # set sample names as index
    train_data.set_index('sample_name', inplace=True)
    test_data.set_index('sample_name', inplace=True)
    return train_data, test_data

def parse_data_new_rf(file_path):
    data = pd.read_csv(file_path, engine='python')
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


def get_hyperparameters_grid():
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


def sample_configurations(hyperparameters_grid, num_of_configurations_to_sample):
    configurations = []
    for i in range(num_of_configurations_to_sample):
        configuration = {}
        for key in hyperparameters_grid:
            configuration[key] = np.random.choice(hyperparameters_grid[key], size=1)[0]
        configurations.append(configuration)
    return configurations


def plot_heat_map(df, number_of_features, output_path, hits, number_of_samples):
    #plt.figure(dpi=3000)
    cm = sns.clustermap(df, cmap="Blues", col_cluster=False, yticklabels=True)
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=150/number_of_samples)
    cm.ax_heatmap.set_title(f"A heat-map of the significance of the top {number_of_features} discriminatory motifs")
    cm.savefig(f"{output_path}/{number_of_features}.svg", format='svg', bbox_inches="tight")
    plt.close()


def plot_error_rate(errors, features, output_path_dir):
    plt.figure(dpi=1000)
    plt.plot(features, errors, '--o')
    plt.xscale('log')
    plt.ylim(-0.02, 1)
    plt.xlabel("features")
    plt.ylabel("error rate")
    plt.savefig(f"{output_path_dir}/error_rate.png")
    plt.close()

def plot_error_rate_new_rf(errors, features, cv_num_of_splits, output_path_dir):
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


def train_models(csv_file_path, num_of_configurations_to_sample, done_path, cv_num_of_splits, use_new_rf, use_tfidf, num_of_iterations, argv):
    logging.info('Parsing data...')
    if use_new_rf:
        X_train, y_train, X_test, y_test, feature_names, sample_names_train, sample_names_test = parse_data_new_rf(csv_file_path)
    else:
        train_data, test_data = parse_data(csv_file_path)
        y = np.array(train_data['label']) # saving the (true) labels
        train_data.drop(['label'], axis=1, inplace=True)
        test_data.drop(['label'], axis=1, inplace=True)
        X = np.array(train_data)  # saving an array of the variables
        max_instances_per_class = np.max(np.unique(y, return_counts=True)[1])

    is_hits_data = 'hits' in csv_file_path

    logging.info('Preparing output path...')
    csv_folder, csv_file_name = os.path.split(csv_file_path)
    csv_file_prefix = os.path.splitext(csv_file_name)[0]  # without extension
    output_path = os.path.join(csv_folder, f'{csv_file_prefix}_model')
    feature_selection_summary_path = f'{output_path}/feature_selection_summary.txt'


    if not use_new_rf:
        for i in range(num_of_iterations):
            output_path_i = os.path.join(output_path, str(i))
            if not os.path.exists(output_path_i):
                logging.info('Creating output path...')
                os.makedirs(output_path_i)

            errors, features = train(X, y, max_instances_per_class, is_hits_data, 
                                    train_data, output_path_i, 3481 + i, use_tfidf)

            plot_error_rate(errors, features, output_path_i)

            with open(feature_selection_summary_path, 'w' if i == 0 else 'a') as f:  # override only on first iteration
                f.write(f'{i}\t{features[-1]}\n')
    else:
        os.makedirs(output_path, exist_ok=True)	
        best_model_path = os.path.join(output_path, f'best_model')	

        # single feature analysis	
        logging.info('Applying single feature analysis...')	
        perfect_feature_names, perfect_feature_indexes = measure_each_feature_accuracy(X_train, y_train, feature_names, output_path)	
        if perfect_feature_names:	
            df = save_model_features(X_train, perfect_feature_indexes, perfect_feature_names, sample_names_train, f'{output_path}/perfect_feature_names')	
            plot_heat_map(df, df.shape[1], output_path, False, df.shape[0])	
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
        hyperparameters_grid = get_hyperparameters_grid()	
        feature_selection_summary_f = open(feature_selection_summary_path, 'w')	
        feature_selection_summary_f.write(f'model_number\tnum_of_features\tfinal_error_rate\n')	
        sampled_configurations = sample_configurations(hyperparameters_grid, num_of_configurations_to_sample)	
        for i, configuration in enumerate(sampled_configurations):	
            model_number = str(i).zfill(len(str(num_of_configurations_to_sample)))	
            output_path_i = os.path.join(output_path, model_number)	
            logging.info(f'Creating output path #{i}...')	
            os.makedirs(output_path_i, exist_ok=True)	
            save_configuration_to_txt_file(configuration, output_path_i)	
            logging.info(f'Configuration #{i} hyper-parameters are:\n{configuration}')	
            rf = RandomForestClassifier(**configuration)	
            errors, features = train_new_rf(rf, X_train, y_train, feature_names, sample_names_train, is_hits_data, output_path_i)	
            plot_error_rate_new_rf(errors, features, cv_num_of_splits, output_path_i)	
            feature_selection_summary_f.write(f'{model_number}\t{features[-1]}\t{errors[-1]}\n')	
            if features[-1] == 1 and errors[-1] == 0:	
                # found the best model (accuracy-wise)	
                # there is no point to continue...	
                break	
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

def save_model_features(X, feature_indexes, feature_names, sample_names, output_path):
    df = pd.DataFrame()
    for i in range(len(feature_names)):
        df[feature_names[i]] = X[:, feature_indexes[i]]  # train_data.iloc[:, :number_of_features]
    df.set_index(sample_names, inplace=True)
    df.to_csv(f"{output_path}.csv")
    return df

def generate_roc_curve(X, y, classifier, number_of_features, output_path):
    ax = plt.gca()
    plot_roc_curve(classifier, X, y, **{'marker': 'o'}, ax=ax)
    plt.savefig(f"{output_path}/roc_curve_{number_of_features}.png", format='png', bbox_inches="tight")
    plt.close()

def train(X, y, max_instances_per_class, hits_data, train_data, output_path, seed, use_tfidf):
    logging.info('Training...')
    rf = RandomForestClassifier(n_estimators=1000, random_state=np.random.seed(seed))  # number of trees
    model = rf.fit(X, y)
    importance = model.feature_importances_
    indexes = np.argsort(importance)[::-1]  # decreasing order of importance
    train_data = train_data.iloc[:, indexes]  # sort features by their importance

    with open(f'{output_path}/feature_importance.txt', 'w') as f:
        importance=importance[indexes]
        features = train_data.columns.tolist()
        for i in range(len(importance)):
            f.write(f'{features[i]}\t{importance[i]}\n')

    # transform the data for better contrast in the visualization
    if hits_data:  # hits data
        train_data = np.log2(train_data+1)  # pseudo counts
    elif use_tfidf: # tfidf data
        train_data = -np.log2(train_data+0.0001)  # avoid 0
    else:  # p-values data
        train_data = -np.log2(train_data)

    number_of_samples, number_of_features = X.shape
    error_rate = previous_error_rate = 1
    error_rates = []
    number_of_features_per_model = []
    while error_rate <= previous_error_rate and number_of_features >= 1:
        logger.info(f'Number of features is {number_of_features}')
        number_of_features_per_model.append(number_of_features)

        # save previous error_rate to make sure the performances do not deteriorate
        previous_error_rate = error_rate

        # compute current model accuracy for each fold of the cross validation
        cv = np.min([max_instances_per_class, 3])
        cv_score = cross_val_score(rf, X, y, cv=cv, n_jobs=-1)

        # current model error_rate rate
        error_rate = 1 - cv_score.mean()
        error_rates.append(error_rate)
        logger.info(f'Error rate is {error_rate}')

        # save current model features to a csv file
        df = train_data.iloc[:, :number_of_features]
        if error_rate <= previous_error_rate:
            df.to_csv(f"{output_path}/Top_{number_of_features}_features.csv")

            plot_heat_map(df, number_of_features, output_path, hits_data, number_of_samples)

            # save the model itself (serialized) for future use
            joblib.dump(model,
                        os.path.join(output_path, f'Top_{number_of_features}_features_model.pkl'))

        # update number of features
        number_of_features //= 2
        if number_of_features < 1:
            continue

        # extract only the (new) half most important features
        X = np.array(train_data.iloc[:, :number_of_features])

        # re-evaluate
        # rf = RandomForestClassifier(n_estimators=100)
        # model = rf.fit(X, y)

        # Sanity check for debugging: predicting the test data
        # change the logging level (line 12) to logging.DEBUG to get the predictions
        # logger.debug(model.predict(test_data.iloc[:, indexes[:number_of_features]]))
    return error_rates, number_of_features_per_model


def train_new_rf(rf, X, y, feature_names, sample_names, hits_data, output_path):
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

        # predictions = model.predict(X)
        # logger.info(f'Full model error rate is {1 - (predictions == y).mean()}')
        # logger.info(f'Current model\'s predictions\n{predictions.tolist()}')

        # compute current model accuracy for each fold of the cross validation
        cv_score = cross_val_score(model, X, y, cv=StratifiedKFold(n_splits=2))

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

        plot_heat_map(df, number_of_features, output_path, hits_data, number_of_samples)

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


def measure_each_feature_accuracy(X_train, y_train, feature_names, output_path):
    feature_to_avg_accuracy = {}
    # df = pd.read_csv(f'{output_path}/Top_149_features.csv', index_col='sample_name')
    rf = RandomForestClassifier()

    for i, feature in enumerate(feature_names):
        # if i % 10 == 0:
        logger.info(f'Checking feature {feature} number {i}')
        # assert df.columns[i] == feature
        cv_score = cross_val_score(rf, X_train[:, i].reshape(-1, 1), y_train, cv=StratifiedKFold(n_splits=2, shuffle=True)).mean()
        if cv_score == 1:
            logger.info('-' * 10 + f'{feature} has 100% accuracy!' + '-' * 10)
        #     print(X_train[:, i].reshape(-1, 1).tolist()[:8] + X_train[:, i].reshape(-1, 1).tolist()[12:])
        #     print(f'min of other class: {min(X_train[:, i].reshape(-1, 1).tolist()[:8] + X_train[:, i].reshape(-1, 1).tolist()[12:])}')
        #     print(X_train[:, i].reshape(-1, 1).tolist()[8:12])
        #     # print(y_train.tolist())
        #     print('Accuracy is 1')
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
    with open(f'{output_path_i}/hyperparameters_configuration.txt', 'w') as f:
        for key in sampled_configuration:
            f.write(f'{key}={sampled_configuration[key]}\n')
    joblib.dump(sampled_configuration, f'{output_path_i}/hyperparameters_configuration.pkl')

if __name__ == '__main__':

    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', type=str, help='A csv file with data matrix to model ')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('num_of_configurations_to_sample', default=1000, type=int, help='How many random configurations of hyperparameters should be sampled?')
    parser.add_argument('--cv_num_of_splits', default=4, help='How folds should be in the cross validation process? (use 0 for leave one out)')
    parser.add_argument('--num_of_iterations', default=10, help='How many should the RF run?')
    parser.add_argument('--new_rf', action='store_true', help='run new random forest version')
    parser.add_argument('--tfidf', action='store_true', help="Are inputs from TF-IDF (avoid log(0))")
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    train_models(args.data_path, args.num_of_configurations_to_sample, args.done_file_path, args.cv_num_of_splits, args.new_rf, args.tfidf, args.num_of_iterations, argv=sys.argv)
