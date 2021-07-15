# running example:
# set exp=Exp26B; python /groups/pupko/orenavr2/gershoni/src/JoinDSSamples.py --domain-name-to-domain-sequence-file /groups/pupko/orenavr2/gershoni/domains_of_interest/HCV-HIV_P3.txt --sample-barcode-to-sample-name-file /groups/pupko/orenavr2/gershoni/Experiments/$exp/data/DS_Samples.txt --lib-name P3 --out-dir ~/gershoni/Experiments/$exp/first_phase_output/

from old.Auxilaries import *
start = time()

logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO) # set level to logging.DEBUG to see debugging comments
logger.info(f'argv: {argv} ; len(argv)= {len(argv)}')

def split_dict_by_interest(domains_of_interest, domain_to_counts, correction):
    of_interest = dict.fromkeys(domains_of_interest, '0') # set to 0 by default
    of_interest_corrected = dict.fromkeys(domains_of_interest, correction) # set to $correction by default
    unrecognized_domains = {}
    for domain in domain_to_counts:
        if domain in domains_of_interest:
            of_interest[domain] = domain_to_counts[domain]
            of_interest_corrected[domain] = str(int(domain_to_counts[domain]) + int(of_interest_corrected[domain]))
            logger.debug('of_interest_corrected[domain] = {}'.format(of_interest_corrected[domain]))
        else:
            unrecognized_domains[domain] = domain_to_counts[domain]
    return of_interest, of_interest_corrected, unrecognized_domains

# if DEFAULT_PARAMS:
#     # logger.error('Usage: python '+argv[0]+' <domain_name_to_domain_sequence_file> <out_dir>; <?correction>')
#     # exit(-1)
#     domain_name_to_domain_sequence_file = '/Users/Oren/Dropbox/Projects/gershoni/domains_of_interest/HCV-HIV_P8.txt'
#     barcode_to_name_file = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/in_test/first_phase/DS_Samples.txt'
#     out_dir = '/Users/Oren/Dropbox/Projects/gershoni/Experiments/out_test/'
#     correction = '50'
#     args.lib_name = os.path.split(domain_name_to_domain_sequence_file)[-1].split('_')[1].split('.')[0]

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--domain-name-to-domain-sequence-file', help='A tab delimited file that maps from domain_name to domain_sequence',
                    type=lambda file_path: str(file_path) if os.path.exists(file_path) else parser.error(
                        f'{file_path} does not exist!'), default='/Users/Oren/Dropbox/Projects/gershoni/domains_of_interest/HCV-HIV_P8.txt')
parser.add_argument('--sample-barcode-to-sample-name-file',
                    help='A tab delimited file that maps from sample sample_barcode to sample name',
                    type=lambda file_path: str(file_path) if os.path.exists(
                        os.path.split(file_path)[0]) else parser.error(
                        f'{os.path.split(file_path)[0]} does not exist!'), default='/Users/Oren/Dropbox/Projects/gershoni/Experiments/in_test/first_phase/DS_Samples.txt')
parser.add_argument('--lib-name', help='Name of the library over which the samples were scanned', default='HIV')
parser.add_argument('--out-dir', help='A path where the output of FilterReads is located (usually called first_phase_output) and to which the output will be written',
                    default='/Users/Oren/Dropbox/Projects/gershoni/Experiments/out_test/')
parser.add_argument('--correction', help='Pseudo counts for normalization',
                    type=lambda x: x if int(x) >= 0 and int(x) == float(x) else parser.error(
                        f'mistakes_allowed should be a positive integer'), default='50')
parser.add_argument('--extended-output', help='Extract additional output files (mainly for sanity checks and debugging)',
                    action='store_true')  # false by default

args = parser.parse_args()

# domain_name_to_domain_sequence_file = args.domain_name_to_domain_sequence_file
# sample_barcode_to_sample_name_file = args.sample_barcode_to_sample_name_file
# out_dir = args.out_dir.rstrip('/')
# correction = args.correction
# lib_name = args.lib_name


logger.info(f'correction={args.correction}')
logger.info(f'lib_name={args.lib_name}')
copy_number_file = os.path.join(args.out_dir, args.lib_name + '.copy_number.csv')
corrected_copy_number_file = os.path.join(args.out_dir, args.lib_name + '.copy_number.corrected_by_'+args.correction+'.csv')
corrected_percents_file = os.path.join(args.out_dir, args.lib_name + '.percent.corrected_by_'+args.correction+'.csv')
recognized_vs_unrecognized_file = os.path.join(args.out_dir, args.lib_name + '.recognized_vs_unrecognized.csv')


domains_of_interest_to_names = load_table_to_dict(args.domain_name_to_domain_sequence_file, backwards=True)
names_to_domains_of_interest = load_table_to_dict(args.domain_name_to_domain_sequence_file)


# list of relevant domains (in the same order they appear in the file!)
domains_of_interest_list = get_column_from_file(args.domain_name_to_domain_sequence_file, 1)
logger.info(f'domains_of_interest_list = {domains_of_interest_list}')

samples_names = get_column_from_file(args.sample_barcode_to_sample_name_file, 1)
logger.info(f'samples_names = {samples_names}')

copy_number_result = ','.join(['sample'] + [domains_of_interest_to_names[domain] for domain in domains_of_interest_list] + ['TOTAL']) + '\n'
corrected_copy_number_result = copy_number_result
corrected_percents_result = ','.join(['sample'] + [domains_of_interest_to_names[domain] for domain in domains_of_interest_list]) + '\n'
recognized_vs_unrecognized_result = ','.join(['recognized_file', 'sum_recognized', 'percent_recognized', 'unrecognized_file', 'sum_unrecognized', 'percent_unrecognized', 'total_reads']) + '\n'

data_is_not_empty = False

logger.info('Aggregating the following files:')
#initializing resources
for sample_name in samples_names:
    sample_path = os.path.join(args.out_dir, sample_name)
    logger.debug(f'sample_path: {sample_path}')
    if not args.lib_name.lower() in sample_name.lower():
        logger.info(f'lib_name: {args.lib_name} NOT in sample_name: {sample_name.lower()}')
    if os.path.isdir(sample_path) and args.lib_name.lower() in sample_name.lower(): # scan only relevant dirs
        for file_name in os.listdir(sample_path):
            logger.debug(f'file_name {file_name}')
            if ALL in file_name:
                counts_file_path = os.path.join(sample_path, file_name)
                domain_to_counts = load_table_to_dict(counts_file_path)
                domain_of_interest_to_counts, domain_of_interest_to_corrected_counts, unrecognized_domains_to_counts =\
                    split_dict_by_interest(set(domains_of_interest_list), domain_to_counts, args.correction)

                logger.info('Aggregating {}'.format(counts_file_path))

                of_interest_path = os.path.join(sample_path, file_name.replace(ALL, DOMAIN_OF_INTEREST))
                with open(of_interest_path, 'w') as f:
                    # domain_of_interest_to_counts keys are exactly domains_of_interest_list but
                    # we need to keep the order in domains_of_interest_list
                    f.write('\n'.join(
                        domains_of_interest_to_names[domain] + '\t' + domain_of_interest_to_counts[domain] for domain in domains_of_interest_list))


                unrecognized_domains_path = os.path.join(sample_path, file_name.replace(ALL, LEFTOVERS))
                sorted_unrecognized_domains = sorted(unrecognized_domains_to_counts, key=lambda x:int(unrecognized_domains_to_counts[x]), reverse=True)
                if args.extended_output:
                    with open(unrecognized_domains_path, 'w') as f:
                        print(unrecognized_domains_path)
                        f.write('\n'.join(
                            domain + '\t' +
                            unrecognized_domains_to_counts[domain] for domain in sorted_unrecognized_domains))

                data_is_not_empty = True # to determine if the files should be written or not

                #copy_number
                total_of_interest = sum(int(domain_of_interest_to_counts[domain]) for domain in domains_of_interest_list)
                copy_number_result += ','.join([sample_name] + [domain_of_interest_to_counts[domain] for domain in domains_of_interest_list] + [str(total_of_interest)]) + '\n'

                #corrected_copy_number
                total_of_interest_corrected = sum(int(domain_of_interest_to_corrected_counts[domain]) for domain in domains_of_interest_list)
                corrected_copy_number_result += ','.join([sample_name] + [str(domain_of_interest_to_corrected_counts[domain]) for domain in domains_of_interest_list] + [str(total_of_interest_corrected)]) + '\n'

                #corrected_percent
                corrected_percents = [float(domain_of_interest_to_corrected_counts[domain]) / total_of_interest_corrected for domain in domains_of_interest_list]
                corrected_percents_result += ','.join([sample_name] + [f'{percent:.7f}' for percent in corrected_percents]) + '\n'

                # find info file path
                for file in os.listdir(sample_path):
                    if INFO in file.lower():
                        info_file_path = os.path.join(sample_path, file)
                        break
                # find total sequence number
                with open(info_file_path) as f:
                    for line in f:
                        if TOTAL_SEQ_AFTER_FILTERS in line:
                            total_seq_after_filters = int(line.rstrip().split()[-1])
                            break

                # recognized_vs_unrecognized
                total_unrecognized_domains = total_seq_after_filters - total_of_interest
                percent_of_interest = total_of_interest / float(total_seq_after_filters) * 100
                percent_unrecognized_domains = total_unrecognized_domains / float(total_seq_after_filters) * 100
                recognized_vs_unrecognized_result += ','.join([of_interest_path, str(total_of_interest), f'{percent_of_interest:.7f}', unrecognized_domains_path, str(total_unrecognized_domains),
                                                               f'{percent_unrecognized_domains:.7f}']) + '\n'

if data_is_not_empty:
    with open(copy_number_file, 'w') as f:
        f.write(copy_number_result)

    with open(corrected_copy_number_file, 'w') as f:
        f.write(corrected_copy_number_result)

    with open(corrected_percents_file, 'w') as f:
        f.write(corrected_percents_result)

    with open(recognized_vs_unrecognized_file, 'w') as f:
        f.write(recognized_vs_unrecognized_result)
else:
    logger.info('No relevant data for this run (no files will be written)')

logger.info('Done aggregating.')
end = time()
logger.info(f'Finished joining samples. Took {measure_time(start, end)}')