import datetime
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *

def run_pipeline(fastq_path, barcode2samplename_path, samplename2biologicalcondition_path, analysis_dir, logs_dir,
                 left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, gz,
                 max_msas_per_sample, max_msas_per_bc,
                 max_number_of_cluster_members_per_sample, max_number_of_cluster_members_per_bc,
                 allowed_gap_frequency, concurrent_cutoffs, meme_split_size, number_of_random_pssms,
                 rank_method, tfidf_method, tfidf_factor, shuffles, use_new_rf, number_of_hyperparams_configuration_to_sample,
                 run_summary_path, error_path, queue, verbose, argv):

    os.makedirs(os.path.split(run_summary_path)[0], exist_ok=True)

    f_run_summary_path = open(run_summary_path, 'w')
    f_run_summary_path.write(' '.join(argv) + '\n\n')
    f_run_summary_path.flush()

    start_time = datetime.datetime.now()

    exp_name = analysis_dir.rstrip('/').split('/')[-2]

    # output folders of the different modules
    first_phase_output_path = os.path.join(analysis_dir, 'reads_filtration')
    second_phase_output_path = os.path.join(analysis_dir, 'motif_inference')
    third_phase_output_path = os.path.join(analysis_dir, 'model_fitting')

    first_phase_done_path = f'{logs_dir}/reads_filtration_done.txt'
    if not os.path.exists(first_phase_done_path):
        os.makedirs(first_phase_output_path, exist_ok=True)
        first_phase_logs_path = os.path.join(logs_dir, 'reads_filtration')
        os.makedirs(first_phase_logs_path, exist_ok=True)

        module_parameters = [fastq_path, first_phase_output_path, first_phase_logs_path,
                             barcode2samplename_path, left_construct, right_construct,
                             max_mismatches_allowed, min_sequencing_quality, first_phase_done_path,
                             '--gz' if gz else '', f'--error_path {error_path}', '-v' if verbose else '']
        cmd = submit_pipeline_step(f'{src_dir}/reads_filtration/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_reads_filtration',
                             queue, verbose)

        wait_for_results('reads_filtration', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='reads_filtration_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping reads filtration. Done file exists at:\n{first_phase_done_path}')

    second_phase_done_path = f'{logs_dir}/motif_inference_done.txt'
    if not os.path.exists(second_phase_done_path):
        os.makedirs(second_phase_output_path, exist_ok=True)
        second_phase_logs_path = os.path.join(logs_dir, 'motif_inference')
        os.makedirs(second_phase_logs_path, exist_ok=True)

        module_parameters = [first_phase_output_path, second_phase_output_path, second_phase_logs_path,
                             samplename2biologicalcondition_path, max_msas_per_sample, max_msas_per_bc,
                             max_number_of_cluster_members_per_sample, max_number_of_cluster_members_per_bc,
                             allowed_gap_frequency, second_phase_done_path,
                             f'--meme_split_size {meme_split_size}',
                             f'--error_path {error_path}', '-v' if verbose else '', f'-q {queue}']
        if concurrent_cutoffs:
            module_parameters.append('--concurrent_cutoffs')
        cmd = submit_pipeline_step(f'{src_dir}/motif_inference/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_motif_inference',
                             queue, verbose)

        wait_for_results('motif_inference', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='motif_inference_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping motif inference. Done file exists at:\n{second_phase_done_path}')

    third_phase_done_path = f'{logs_dir}/model_fitting_done.txt'
    if not os.path.exists(third_phase_done_path):
        os.makedirs(third_phase_output_path, exist_ok=True)
        third_phase_logs_path = os.path.join(logs_dir, 'model_fitting')
        os.makedirs(third_phase_logs_path, exist_ok=True)

        module_parameters = [first_phase_output_path, second_phase_output_path, third_phase_output_path,
                             third_phase_logs_path, samplename2biologicalcondition_path, number_of_random_pssms, number_of_hyperparams_configuration_to_sample,
                             third_phase_done_path, f'--rank_method {rank_method}', f'--error_path {error_path}', '--new_rf' if use_new_rf else ''
                             '-v' if verbose else '', f'-q {queue}']
        if rank_method == 'tfidf':
            if tfidf_method:
                module_parameters += ['--tfidf_method', tfidf_method]
            if tfidf_factor:
                module_parameters += ['--tfidf_factor', str(tfidf_factor)]
        elif rank_method == 'shuffles':
            if shuffles:
                module_parameters += ['--shuffles', shuffles]
        cmd = submit_pipeline_step(f'{src_dir}/model_fitting/module_wraper.py',
                             [module_parameters],
                             logs_dir, f'{exp_name}_model_fitting',
                             queue, verbose)

        wait_for_results('model_fitting', logs_dir, num_of_expected_results=1, example_cmd=cmd,
                         error_file_path=error_path, suffix='model_fitting_done.txt')
    else:
        logger.info(f'{datetime.datetime.now()}: skipping model fitting. Done file exists {third_phase_done_path}')

    end_time = datetime.datetime.now()
    f_run_summary_path.write(f'Total running time: {str(end_time-start_time)[:-3]}')
    f_run_summary_path.close()

    logger.info(f'Started running at {start_time}')
    logger.info(f'Done running at {end_time}')
    logger.info(f'Total running time: {str(end_time-start_time)[:-3]}')
    logger.info('Bye!')



if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}', flush=True)

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('barcode2samplename_path', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('samplename2biologicalcondition_path', type=str, help='A path to the sample name to biological condition file')
    parser.add_argument('analysis_dir', type=str, help='analysis folder')
    parser.add_argument('logs_dir', type=str, help='logs folder')

    # optional parameters for the filtration step
    parser.add_argument('--left_construct', default='CAACGTGGC', help='The (constant) sequence from the LEFT of the random sequence') # in exp12: "CAACGTGGC"
    parser.add_argument('--right_construct', default='GCCT', help='The (constant) sequence from the RIGHT of the random sequence') # in exp12: "GCCT"
    parser.add_argument('--max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constant sequences')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')

    # optional parameters for the motif inference
    parser.add_argument('--max_msas_per_sample', default=100, type=int,
                        help='For each sample, align only the biggest $max_msas_per_sample')
    parser.add_argument('--max_msas_per_bc', default=400, type=int,
                        help='For each biological condition, align only the biggest $max_msas_per_bc')
    parser.add_argument('--max_number_of_cluster_members_per_sample', default=1000, type=int,
                        help='How many members (at most) should be taken to each cluster')
    parser.add_argument('--max_number_of_cluster_members_per_bc', default=100, type=int,
                        help='How many members (at most) should be taken to each cluster')
    parser.add_argument('--allowed_gap_frequency', default=0.9,
                        help='Maximal gap frequency allowed in msa (higher frequency columns are removed)',
                        type=lambda x: float(x) if 0 < float(x) < 1
                                                else parser.error(f'The threshold of the maximal gap frequency allowed per column should be between 0 to 1'))
    parser.add_argument('--concurrent_cutoffs', action='store_true',
                        help='Use new method which splits meme before cutoffs and runs cutoffs concurrently')
    parser.add_argument('--meme_split_size', type=int, default=1, # TODO default of 1, 5 or 10?
                        help='Split size, how many meme per files for calculations')

    # optional parameters for the modelling step
    parser.add_argument('--number_of_random_pssms', default=100, type=int, help='Number of pssm permutations')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles'], default='pval', help='Motifs ranking method')
    parser.add_argument('--tfidf_method', choices=['boolean', 'terms', 'log', 'augmented'], default='boolean', help='TF-IDF method')
    parser.add_argument('--tfidf_factor', type=float, default=0.5, help='TF-IDF augmented method factor (0-1)')
    parser.add_argument('--shuffles', default=5, type=int, help='Number of controlled shuffles permutations')
    parser.add_argument('--number_of_hyperparams_configuration_to_sample', default=1000, type=int,
                        help='How many random configurations of hyperparameters should be sampled for choosing the RF model?')

    # general optional parameters
    parser.add_argument('--run_summary_path', type=str,
                        help='a file in which the running configuration and timing will be written to')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
    parser.add_argument('--new_rf', action='store_true', help='run new random forest version')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    run_summary_path = args.error_path if args.error_path else os.path.join(args.analysis_dir, 'run_summary_path.txt')
    error_path = args.error_path if args.error_path else os.path.join(args.logs_dir, 'error.txt')

    concurrent_cutoffs = True if args.concurrent_cutoffs else False

    run_pipeline(args.fastq_path, args.barcode2samplename_path, args.samplename2biologicalcondition_path,
                 args.analysis_dir.rstrip('/'), args.logs_dir.rstrip('/'),
                 args.left_construct, args.right_construct, args.max_mismatches_allowed, args.min_sequencing_quality, True if args.gz else False,
                 args.max_msas_per_sample, args.max_msas_per_bc,
                 args.max_number_of_cluster_members_per_sample, args.max_number_of_cluster_members_per_bc,
                 args.allowed_gap_frequency, concurrent_cutoffs, args.meme_split_size, args.number_of_random_pssms,
                 args.rank_method, args.tfidf_method, args.tfidf_factor, args.shuffles, args.new_rf, args.number_of_hyperparams_configuration_to_sample,
                 run_summary_path, error_path, args.queue, True if args.verbose else False, sys.argv)

