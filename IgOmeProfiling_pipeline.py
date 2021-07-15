import datetime
from genericpath import exists
import os
import sys
import json 
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

from auxiliaries.pipeline_auxiliaries import *
from auxiliaries.validation_files import is_input_files_valid
from auxiliaries.stop_machine_aws import stop_machines

schema = {
    'reads': schema_reads,
    'inference': schema_inference,
    'model': schema_cross_exp
}

def validate_input_file(config, phase, logger):
    f = open(config)
    config_dict = json.load(f)
    samplename2biologicalcondition_paths = []
    barcode2samplename_paths = []
    is_valid = is_valid_json_structure(config, config_dict, schema[phase], logger)
    if not is_valid:
        logger.info(f'Json file {config} is not valid')
        return False
    runs = config_dict['runs']
    for run in runs:
        if phase == 'reads':
            barcode2samplename_paths.append(runs[run]['barcode2sample'])
        if phase == 'inference':
            samplename2biologicalcondition_paths.append(runs[run]['sample2bc'])
        if phase == 'model':
            samplename2biologicalcondition_paths.append(runs[run]['sample2bc'].values())
    for path in barcode2samplename_paths:
        is_valid = is_input_files_valid(samplename2biologicalcondition_path='', barcode2samplename_path=path, logger=logger)
        if not is_valid:
            return False
    samplename2biologicalcondition_paths = list(set(samplename2biologicalcondition_paths))
    for path in samplename2biologicalcondition_paths:
        is_valid = is_input_files_valid(samplename2biologicalcondition_path=path, barcode2samplename_path='', logger=logger)
        if not is_valid:
            return False        
    return True


def is_all_files_valid(multi_exp_config_reads, multi_exp_config_inference, cross_experiments_config, samplename2biologicalcondition_path, barcode2samplename_path, logger):
    is_files_valid = True
    is_files_valid_reads = True
    is_files_valid_inference = True
    is_files_valid_model = True

    if os.path.exists(samplename2biologicalcondition_path) and os.path.exists(barcode2samplename_path):
        is_files_valid = is_input_files_valid(samplename2biologicalcondition_path=samplename2biologicalcondition_path, barcode2samplename_path=barcode2samplename_path, logger=logger)
    else:
        if multi_exp_config_reads:
            is_files_valid_reads = validate_input_file(multi_exp_config_reads, 'reads', logger)
        if multi_exp_config_inference:
            is_files_valid_inference = validate_input_file(multi_exp_config_inference, 'inference', logger)
        if cross_experiments_config:
            is_files_valid_model = validate_input_file(cross_experiments_config, 'model', logger)

    if not is_files_valid or not is_files_valid_reads or not is_files_valid_inference or not is_files_valid_model:
        return False
    return True


def run_pipeline(fastq_path, barcode2samplename_path, samplename2biologicalcondition_path, analysis_dir, logs_dir,
                 left_construct, right_construct, max_mismatches_allowed, min_sequencing_quality, minimal_length_required, gz, rpm, multi_exp_config_reads,
                 max_msas_per_sample, max_msas_per_bc, max_number_of_cluster_members_per_sample, max_number_of_cluster_members_per_bc,
                 allowed_gap_frequency, threshold, word_length, discard, cluster_algorithm_mode, concurrent_cutoffs, meme_split_size, use_mapitope, aln_cutoff,
                 pcc_cutoff, skip_sample_merge_meme, minimal_number_of_columns_required_create_meme, prefix_length_in_clstr, multi_exp_config_inference, cutoff_random_peptitdes_percentile,
                 stop_before_random_forest, is_run_random_forest_per_bc_sequentially, number_of_random_pssms, number_parallel_random_forest, min_value_error_random_forest,
                 rank_method, tfidf_method, tfidf_factor, shuffles, shuffles_percent, shuffles_digits,
                 num_of_random_configurations_to_sample, cv_num_of_splits, seed_random_forest, random_forest_seed_configurations,
                 stop_machines_flag, type_machines_to_stop, name_machines_to_stop, cross_experiments_config, no_rpm_factor,
                 run_summary_path, error_path, queue, verbose, argv):

    files_are_valid = True
    if multi_exp_config_reads or multi_exp_config_inference or cross_experiments_config:
        files_are_valid = is_all_files_valid(multi_exp_config_reads, multi_exp_config_inference, cross_experiments_config,
                                             samplename2biologicalcondition_path, barcode2samplename_path, logger)
    else:    
    # check the validation of files barcode2samplename_path and samplename2biologicalcondition_path
        files_are_valid = is_input_files_valid(samplename2biologicalcondition_path=samplename2biologicalcondition_path, barcode2samplename_path=barcode2samplename_path, logger=logger)
    if not files_are_valid:
        return

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
                            max_mismatches_allowed, min_sequencing_quality, first_phase_done_path, minimal_length_required, '--check_files_valid' if not files_are_valid else '',
                            f'--multi_exp_config_reads {multi_exp_config_reads}' if multi_exp_config_reads else '',
                            '--rpm' if rpm else '', '--gz' if gz else '', f'--error_path {error_path}', '-v' if verbose else '', '-m' if use_mapitope else '']        
        
        cmd = submit_pipeline_step(f'{src_dir}/reads_filtration/module_wraper.py',[module_parameters],
                             logs_dir, f'{exp_name}_reads_filtration', queue, verbose)

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
                             allowed_gap_frequency, second_phase_done_path, '--check_files_valid' if not files_are_valid else '', 
                             f'--minimal_number_of_columns_required_create_meme {minimal_number_of_columns_required_create_meme}',
                             f'--prefix_length_in_clstr {prefix_length_in_clstr}', f'--aln_cutoff {aln_cutoff}', f'--pcc_cutoff {pcc_cutoff}',
                             f'--threshold {threshold}', f'--word_length {word_length}', f'--discard {discard}', f'--cluster_algorithm_mode {cluster_algorithm_mode}',
                             f'--multi_exp_config_inference {multi_exp_config_inference}' if multi_exp_config_inference else '',
                             f'--meme_split_size {meme_split_size}', f'--skip_sample_merge_meme {skip_sample_merge_meme}',
                             f'--cutoff_random_peptitdes_percentile {cutoff_random_peptitdes_percentile}',
                             f'--error_path {error_path}', '-v' if verbose else '', f'-q {queue}','-m' if use_mapitope else '']       
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
                             third_phase_logs_path, samplename2biologicalcondition_path, number_of_random_pssms, third_phase_done_path,
                             '--stop_before_random_forest' if stop_before_random_forest else '', f'--cross_experiments_config {cross_experiments_config}' if cross_experiments_config else '', 
                             f'--num_of_random_configurations_to_sample {num_of_random_configurations_to_sample}', '--check_files_valid' if not files_are_valid else '',
                             f'--number_parallel_random_forest {number_parallel_random_forest}', f'--min_value_error_random_forest {min_value_error_random_forest}',
                             f'--shuffles_percent {shuffles_percent}', f'--shuffles_digits {shuffles_digits}',
                             f'--cv_num_of_splits {cv_num_of_splits}', f'--seed_random_forest {seed_random_forest}',
                             f'--random_forest_seed_configurations {random_forest_seed_configurations}', f'--rank_method {rank_method}', 
                             '--is_run_random_forest_per_bc_sequentially' if is_run_random_forest_per_bc_sequentially else '',
                             '' if no_rpm_factor else '--no_rpm_factor',
                             f'--error_path {error_path}', '-v' if verbose else '', f'-q {queue}','-m' if use_mapitope else '']        
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

    if stop_machines_flag:
        stop_machines(type_machines_to_stop, name_machines_to_stop, logger)

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
    parser.add_argument('--minimal_length_required', default=3, type=int,
                        help='Shorter peptides will be discarded')                             
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')
    parser.add_argument('--rpm', action='store_true', help='Normalize counts to "reads per million" (sequence proportion x 1,000,000)')
    parser.add_argument('--multi_exp_config_reads', type=str, help='Configuration file for reads phase to run multi expirements')

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
    parser.add_argument('--threshold', default='0.5', help='Minimal sequence similarity threshold required',
                        type=lambda x: float(x) if 0.4 <= float(x) <= 1
                                                else parser.error(f'CD-hit allows thresholds between 0.4 to 1'))
    parser.add_argument('--word_length', default='2', choices=['2', '3', '4', '5'],
                        help='A heuristic of CD-hit. Choose of word size:\n5 for similarity thresholds 0.7 ~ 1.0\n4 for similarity thresholds 0.6 ~ 0.7\n3 for similarity thresholds 0.5 ~ 0.6\n2 for similarity thresholds 0.4 ~ 0.5')
    parser.add_argument('--discard', default='1', help='Include only sequences longer than <$discard> for the analysis. (CD-hit uses only sequences that are longer than 10 amino acids. When the analysis includes shorter sequences, this threshold should be lowered. Thus, it is set here to 1 by default.)')
    parser.add_argument('--cluster_algorithm_mode', default='0', type=str, help='0 - clustered to the first cluster that meet the threshold (fast). 1 - clustered to the most similar cluster (slow)')
    parser.add_argument('--concurrent_cutoffs', action='store_true',
                        help='Use new method which splits meme before cutoffs and runs cutoffs concurrently')
    parser.add_argument('--meme_split_size', type=int, default=1, # TODO default of 1, 5 or 10?
                        help='Split size, how many meme per files for calculations')
    parser.add_argument('-m', '--mapitope', action='store_true', help='use mapitope encoding')
    parser.add_argument('--aln_cutoff', default='20', help='The cutoff for pairwise alignment score to unite motifs of BC') 
    parser.add_argument('--pcc_cutoff', default='0.6', help='Minimal PCC R to unite motifs of BC')
    parser.add_argument('--skip_sample_merge_meme', default='a_weird_str_that_shouldnt_be_a_sample_name_by_any_chance',
                        help='A sample name that should be skipped in merge meme files, e.g., for testing purposes. More than one sample '
                             'name should be separated by commas but no spaces. '
                             'For example: 17b_05,17b_05_test,another_one')
    parser.add_argument('--minimal_number_of_columns_required_create_meme', default=1, type=int,
                        help='MSAs with less than the number of required columns will be skipped')
    parser.add_argument('--prefix_length_in_clstr', default=20, type=int,
                        help='How long should be the prefix that is taken from the clstr file (cd-hit max prefix is 20)')
    parser.add_argument('--multi_exp_config_inference', type=str, help='Configuration file for inference motifs phase to run multi expirements')
    parser.add_argument('--cutoff_random_peptitdes_percentile', type=float, default=0.05, help='Calculate cutoff (hit threshold) from random peptides\' top percentile score')

    # optional parameters for the modelling step
    parser.add_argument('--cross_experiments_config', type=str, help='Configuration file to run cross expiremets at model fitting phase')
    parser.add_argument('--stop_before_random_forest', action='store_true', help='A boolean flag for mark if we need to run the random forest')
    parser.add_argument('--number_of_random_pssms', default=100, type=int, help='Number of pssm permutations')
    parser.add_argument('--number_parallel_random_forest', default=20, type=int, help='How many random forest configurations to run in parallel')
    parser.add_argument('--min_value_error_random_forest', default=0, type=float, help='A random forest model error value for convergence allowing to stop early')
    parser.add_argument('--rank_method', choices=['pval', 'tfidf', 'shuffles'], default='shuffles', help='Motifs ranking method')
    parser.add_argument('--tfidf_method', choices=['boolean', 'terms', 'log', 'augmented'], default='boolean', help='TF-IDF method')
    parser.add_argument('--tfidf_factor', type=float, default=0.5, help='TF-IDF augmented method factor (0-1)')
    parser.add_argument('--shuffles', default=5, type=int, help='Number of controlled shuffles permutations')
    parser.add_argument('--shuffles_percent', default=0.2, type=float, help='Percent from shuffle with greatest number of hits (0-1)')
    parser.add_argument('--shuffles_digits', default=2, type=int, help='Number of digits after the point to print in scanning files.')
    parser.add_argument('--num_of_random_configurations_to_sample', default=100, type=int, help='How many random configurations of hyperparameters should be sampled?')
    parser.add_argument('--cv_num_of_splits', type=int, default=2, help='How folds should be in the cross validation process? (use 0 for leave one out)')
    parser.add_argument('--seed_random_forest', type=int, default=42, help='Seed number for reconstructing experiments')
    parser.add_argument('--random_forest_seed_configurations', default=123 , type=int, help='Random seed value for generating random forest configurations')
    parser.add_argument('--stop_machines', action='store_true', help='Turn off the machines in AWS at the end of the running')
    parser.add_argument('--type_machines_to_stop', default='', type=str, help='Type of machines to stop, separated by comma. Empty value means all machines. Example: t2.2xlarge,m5a.24xlarge ')
    parser.add_argument('--name_machines_to_stop', default='', type=str, help='Names (patterns) of machines to stop, separated by comma. Empty value means all machines. Example: worker*')
    parser.add_argument('--is_run_random_forest_per_bc_sequentially', action='store_true', help='Set the flag to true when number of cores is less than number of BC X 2 (hit and value), otherwise it will run all the BC  parallel (on the same time)')
    parser.add_argument('--no_rpm_factor', action='store_false', help='Disable multiplication hits by factor rpm for normalization')

    # general optional parameters
    parser.add_argument('--run_summary_path', type=str,
                        help='a file in which the running configuration and timing will be written to')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('-q', '--queue', default='pupkoweb', type=str, help='a queue to which the jobs will be submitted')
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
                 args.left_construct, args.right_construct, args.max_mismatches_allowed, args.min_sequencing_quality, args.minimal_length_required,
                 args.gz, args.rpm, args.multi_exp_config_reads,
                 args.max_msas_per_sample, args.max_msas_per_bc, args.max_number_of_cluster_members_per_sample, args.max_number_of_cluster_members_per_bc,
                 args.allowed_gap_frequency, args.threshold, args.word_length, args.discard, args.cluster_algorithm_mode, concurrent_cutoffs, args.meme_split_size, 
                 args.mapitope, args.aln_cutoff, args.pcc_cutoff, args.skip_sample_merge_meme, args.minimal_number_of_columns_required_create_meme,
                 args.prefix_length_in_clstr, args.multi_exp_config_inference, args.cutoff_random_peptitdes_percentile,
                 args.stop_before_random_forest, args.is_run_random_forest_per_bc_sequentially, args.number_of_random_pssms,
                 args.number_parallel_random_forest, args.min_value_error_random_forest,
                 args.rank_method, args.tfidf_method, args.tfidf_factor, args.shuffles, args.shuffles_percent, args.shuffles_digits,
                 args.num_of_random_configurations_to_sample, args.cv_num_of_splits, args.seed_random_forest, args.random_forest_seed_configurations,
                 args.stop_machines, args.type_machines_to_stop, args.name_machines_to_stop, args.cross_experiments_config, args.no_rpm_factor,
                 run_summary_path, error_path, args.queue, args.verbose, sys.argv)
