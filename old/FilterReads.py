#in order to submit this program as a job to the queue, create a copy of the file:
# /groups/pupko/orenavr2/gershoni/Experiments/Exp26A/FilterReads.sge and edit it (change the experiment number) accordingly
    # set exp=Exp26A; python /groups/pupko/orenavr2/gershoni/src/FilterReads.py --fastq-path /groups/pupko/orenavr2/gershoni/Experiments/$exp/data/$exp.fastq.gz --sample-barcode-to-sample-name-file /groups/pupko/orenavr2/gershoni/Experiments/$exp/data/DS_Samples.txt --out-dir /groups/pupko/orenavr2/gershoni/Experiments/$exp/first_phase_output/

# The structure is:
# ===============================================
# Sample barcode + domain upstream sequence + domain barcode + domain downstream sequence
# In the classic Domain Scan experiments it was as follows:
# Sample barcode (5 bp) - no mistakes allowed
# fth1 annealed Site (16 bp) - with mistakes allowed
# Domain barcode (12 bp)
# fth1 annealed - anti sense (17 bp) - no mistakes allowed


from get_args_from_user import *

from old.Auxilaries import *

# sample_barcode_length is defined by the user...
args.barcode_upstream_sequence_length = len(args.barcode_upstream_sequence)  # upstream sequence  is defined by the user...
args.domain_barcode_length = 12
args.barcode_downstream_sequence_length = len(args.barcode_downstream_sequence)  # downstream sequence  is defined by the user...

start = time()

logger.info(f'running command:\npython {" ".join(argv)}\n')

logger = logging.getLogger('main')
logging.basicConfig(level=logging.INFO) # set level to logging.DEBUG to see debugging comments


def filter_fastq(fastq_file, counters, sequence_without_barcode_to_counts):
    local_num_of_seqs_read = 0
    local_num_of_no_barcode = 0
    for record in ParseFastQ(fastq_file):
        # every read is represented by a text chunk of 4 lines.
        # The sequence is in item of the record (4 lines each) e.g.:
        # '@D00257:291:CB65BANXX:4:2310:3171:2243 1:N:0:1\n
        # GGATCAAGTAGGGGATCCAGGTGGTGCGCGACCTCTAGAGCCGACCGCGATAGATCGGAAG\n
        # +\n
        # CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n'

        local_num_of_seqs_read += 1
        if local_num_of_seqs_read % 100000 == 0:
            logger.info(f'num_of_seqs_read is {local_num_of_seqs_read}')

        # next 4-lines chunk
        seq_header, sequence, quality_header, quality_str = record

        sample_barcode = sequence[:args.sample_barcode_length]

        if sample_barcode not in sample_barcode_to_sample_name:
            local_num_of_no_barcode += 1
            sequence_without_barcode_to_counts[sequence] = sequence_without_barcode_to_counts.get(sequence, 0) + 1
            continue  # No sample_barcode -> Skip record...

        # Reached here -> sample_barcode exists.
        putative_barcode_upstream_sequence = sequence[args.sample_barcode_length:args.sample_barcode_length + args.barcode_upstream_sequence_length]
        counters[sample_barcode][TOTAL_SEQ] += 1
        barcode_upstream_sequence_regex = '(' + args.barcode_upstream_sequence + '){s<=' + args.mistakes_allowed + '}'
        if not regex.match(barcode_upstream_sequence_regex, putative_barcode_upstream_sequence):
            counters[sample_barcode][TOO_MANY_MISTAKES] += 1
            N_count = putative_barcode_upstream_sequence.count('N')  # reflects somehow the sequencing quality
            if N_count:  # unrecognized nucleotides
                counters[sample_barcode][WRONG_SEQ_WITH_N] += 1
            if args.extended_output:
                files_handlers[sample_barcode][FILTERED].write('\t'.join([TOO_MANY_MISTAKES, sequence, 'N:', str(N_count), 'START MISTAKES',
                                                                      str(diff(args.barcode_upstream_sequence, putative_barcode_upstream_sequence))]) + '\n')
            continue  # too many mismatches

        # Reached here -> The upstream domain barcode sequence is as expected or with no more than mistakes allowed
        putative_barcode_downstream_sequence = sequence[barcode_upstream_sequence_start_position:barcode_upstream_sequence_start_position + args.barcode_downstream_sequence_length]
        if not args.barcode_downstream_sequence == putative_barcode_downstream_sequence:
            counters[sample_barcode][SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE] += 1
            N_count = putative_barcode_downstream_sequence.count('N')  # reflects the sequencing quality
            if N_count:
                counters[sample_barcode][WRONG_SEQ_WITH_N] += 1
            if args.extended_output:
                files_handlers[sample_barcode][FILTERED].write('\t'.join(['NO_FTH1_Anealed_AntiSense', sequence, 'N:', str(N_count)]) + '\n')
            continue  # Not a PERFECT match

        # Reached here -> The downstream domain barcode sequence matches exactly
        # all tests are OK!!!
        counters[sample_barcode][TOTAL_SEQS_OK] += 1
        domain_barcode = sequence[domain_barcode_start_position:barcode_upstream_sequence_start_position]
        # files_handlers[sample_barcode][FS].write('_'.join(['>Seq', str(counters[sample_barcode][TOTAL_SEQ]), 'LIB_DS']) + '\n')
        files_handlers[sample_barcode][FS].write(f'{domain_barcode}\n')
        if args.extended_output:
            files_handlers[sample_barcode][FASTQ].write(('\n'.join(record)+'\n').encode('utf-8'))
    return local_num_of_seqs_read, local_num_of_no_barcode


domain_barcode_start_position = args.sample_barcode_length + args.barcode_upstream_sequence_length
barcode_upstream_sequence_start_position = domain_barcode_start_position + args.domain_barcode_length


done_file = os.path.join(args.out_dir,'done.txt')
if os.path.exists(done_file):
    logger.info('Results for Exp already exist at:\n' + args.out_dir)
    exit()

counter_types = [TOO_MANY_MISTAKES, SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE, WRONG_SEQ_WITH_N, TOTAL_SEQ, TOTAL_SEQS_OK]

sample_barcode_to_sample_name = load_table_to_dict(args.sample_barcode_to_sample_name_file)

# samples sample_barcodes in the SAME ORDER at the DS_samples file
sample_barcodes = get_column_from_file(args.sample_barcode_to_sample_name_file, 0)

logger.info(sample_barcode_to_sample_name)

files_paths = {}
files_handlers = {}

if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)

num_of_opened_files = 0

logger.info(f'num of sample_barcodes is {len(sample_barcodes)}. {5 if args.extended_output else 3} files will be opened for each sample_barcode ({3 if args.extended_output else 1} now + 2 at the end).')
#initializing resources
for sample_barcode in sample_barcodes:

    name = sample_barcode_to_sample_name.get(sample_barcode, sample_barcode)
    logger.info(f'Barcode: {sample_barcode}')

    # create sub directories for each sample_barcode
    sub_dir = os.path.join(args.out_dir, name)
    if not os.path.exists(sub_dir):
        os.makedirs(sub_dir)

    # initializing file names and handlers for each out_file_type
    files_paths[sample_barcode] = {}
    files_handlers[sample_barcode] = {}

    #file for domain sample_barcodes that passed the QA
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', args.mistakes_allowed, FS])) + args.output_files_suffix
    files_paths[sample_barcode][FS] = out_file
    files_handlers[sample_barcode][FS] = open(out_file, 'w')

    #info file
    #NO need to open it now.
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', args.mistakes_allowed, INFO])) + args.output_files_suffix
    files_paths[sample_barcode][INFO] = out_file
    #files_handlers[sample_barcode][INFO] = open(out_file, 'w')

    #counts file
    #NO need to open it now.
    out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', args.mistakes_allowed, ALL, COUNTS])) + args.output_files_suffix
    files_paths[sample_barcode][COUNTS] = out_file

    if args.extended_output:
        #file for filtered sequences
        out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', args.mistakes_allowed, FILTERED])) + args.output_files_suffix
        files_paths[sample_barcode][FILTERED] = out_file
        files_handlers[sample_barcode][FILTERED] = open(out_file, 'w')

        #file for filtered + domain sample_barcodes
        out_file = os.path.join(sub_dir, '.'.join(['MistakeAllowed', args.mistakes_allowed, FASTQ])) + '.gz'
        files_paths[sample_barcode][FASTQ] = out_file
        files_handlers[sample_barcode][FASTQ] = gzip.open(out_file, 'wb') # these are very (!) large files. Better be kept zipped..


    num_of_opened_files += 3 if args.extended_output else 1
    logger.info('{} were opened.'.format(num_of_opened_files))

    if len(files_handlers) != num_of_opened_files/(3 if args.extended_output else 1):
        logger.error('\n{"#"*50}\nNot enough files were opened. This could be due to duplicated sample_barcodes (i.e., non unique sample_barcodes)\n{"#"*50}')
        exit()

# counters
num_of_no_barcode = 0
num_of_seqs_read = 0
counters = {}
summary_counters = {}

# per sample_barcode counters
for sample_barcode in sample_barcodes:
    counters[sample_barcode] = {}
    for counter in counter_types:
        counters[sample_barcode][counter] = 0

for counter in counter_types:
    summary_counters[counter] = 0

logger.info('Done initializing parameters! Starting to parse the fastq file...')


logger.info('Start parsing the deep sequencing file(s)...')
sequence_without_barcode_to_counts = {}
fastq_files = []
if os.path.isdir(args.fastq_path):
    fastq_files.extend(os.path.join(args.fastq_path, file_name) for file_name in os.listdir(args.fastq_path) if 'fastq' in file_name)
else: #fastq_path is a file
    fastq_files.append(args.fastq_path)
for fastq_file in fastq_files:
    logger.info('\n' + '#'*50 + '\nStart parsing {}\n'.format(fastq_file) + '#'*50)
    local_values = filter_fastq(fastq_file, counters, sequence_without_barcode_to_counts)
    num_of_seqs_read += local_values[0]
    num_of_no_barcode += local_values[1]
logger.info('Done parsing fastq file!')

# Write sequences without barcodes (only when extended_output param is ON)
if args.extended_output:
    no_barcodes_path = os.path.join(args.out_dir, 'no_barcode.txt')
    logger.info(f'Writing sequences without barcodes to: {no_barcodes_path}')
    sorted_sequences_by_counts = sorted(sequence_without_barcode_to_counts, key=sequence_without_barcode_to_counts.get,
                                        reverse=True)
    text_to_no_barcodes_file = ''
    for sequence in sorted_sequences_by_counts:
        text_to_no_barcodes_file += sequence + '\t' + str(sequence_without_barcode_to_counts[sequence]) + '\n'
    with open(no_barcodes_path, 'w') as f:
        f.write(text_to_no_barcodes_file)

logger.info('Writing a summary file for each sample')
for sample_barcode in sample_barcodes:
    text_to_info_file = f'{" ".join(argv)}\n'
    text_to_info_file += f'IN: {args.fastq_path}\n'
    text_to_info_file += f'OUT: {args.out_dir}\n'
    text_to_info_file += f'MISTAKES ALLOWED: {args.mistakes_allowed}\n'
    text_to_info_file += f'BARCODE: {sample_barcode}\n'
    text_to_info_file += f'Total seq: {counters[sample_barcode][TOTAL_SEQ]}\n'
    text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
    text_to_info_file += f'Sequences with more than {args.mistakes_allowed} mistakes in fth1 annealed site sequence region: ' + str(counters[sample_barcode][TOO_MANY_MISTAKES]) + '\n'
    text_to_info_file += f'Wrong seq with N: {counters[sample_barcode][WRONG_SEQ_WITH_N]}\n'
    text_to_info_file += f'Seq with errors in the fth1 annealed anti sense construct: {counters[sample_barcode][SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE]}\n'
    text_to_info_file += f'{TOTAL_SEQ_AFTER_FILTERS}: {counters[sample_barcode][TOTAL_SEQS_OK]}\n'
    with open(files_paths[sample_barcode][INFO], 'w') as f:
        f.write(text_to_info_file)

    # sum up counters of all samples for the general summary file
    for counter in counter_types:
        summary_counters[counter] += counters[sample_barcode][counter]

logger.info('Writing a general summary info file...')
text_to_info_file = f'{" ".join(argv)}\n'
text_to_info_file += f'IN: {args.fastq_path}\n'
text_to_info_file += f'OUT: {args.out_dir}\n'
text_to_info_file += f'MISTAKES ALLOWED: {args.mistakes_allowed}\n'
text_to_info_file += f'BARCODES: {",".join(sample_barcode_to_sample_name)}\n'
text_to_info_file += '\n\n=============================================\n'
text_to_info_file += f'TOTAL ANALYZED SEQ: {num_of_seqs_read}\n'
text_to_info_file += f'SKIPPED SEQ WITH NO BARCODE: {num_of_no_barcode}\n'
for sample_barcode in sample_barcodes:
    text_to_info_file += f'{sample_barcode}: {counters[sample_barcode][TOTAL_SEQ]}\n'
text_to_info_file += f'\nTotal seq with sample_barcodes: {summary_counters[TOTAL_SEQ]}\n'
text_to_info_file += '\n\nWRONG SEQ STATISTICS\n=============================================\n'
text_to_info_file += f'Sequences with more than {args.mistakes_allowed} mistakes in fth1 annealed site sequence region: {summary_counters[TOO_MANY_MISTAKES]}\n'
text_to_info_file += f'Wrong seq with N: {summary_counters[WRONG_SEQ_WITH_N]}\n'
text_to_info_file += f'Seq with no stop codon construct: {summary_counters[SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE]}\n'
text_to_info_file += '\n\nTotal sequences after filters - by sample_barcode:\n=============================================\n'
for sample_barcode in sample_barcodes:
    text_to_info_file += f'{sample_barcode}: {counters[sample_barcode][TOTAL_SEQS_OK]}\n'
text_to_info_file += f'{TOTAL_SEQ_AFTER_FILTERS}: {summary_counters[TOTAL_SEQS_OK]}\n'
with open(os.path.join(args.out_dir, INFO+args.output_files_suffix), 'w') as f:
    f.write(text_to_info_file)

logger.info('Freeing up resources...')
for sample_barcode in files_handlers:
    for handler in files_handlers[sample_barcode]:
        files_handlers[sample_barcode][handler].close()

#frequency_counters_per_sample = {}
logger.info('Writing unique files...')
for sample_barcode in sample_barcodes:
    logger.info('Barcode: ' + sample_barcode)
    #frequency_counters_per_sample[sample_barcode] = {}
    with open(files_paths[sample_barcode][FS]) as f:
        txt = f.read()
    with open(files_paths[sample_barcode][COUNTS], 'w') as f:
        counts = ''
        frequency_counter = Counter(txt.split())
        # with open(files_paths[sample_barcode][COUNTS][:-4]+'2'+files_paths[sample_barcode][COUNTS][-4:]) as f2:
        #     txt = f2.read()
        #     lines = [x.split('_')[-1].strip().split()[::-1] for x in txt.split('>')[1:]]
        #     d = Counter(dict([(x, int(y)) for x,y in lines]))
        #     print d.most_common()
        #     print frequency_counter.most_common()
        #     assert d==frequency_counter
        for seq, count in frequency_counter.most_common():
            counts += f'{seq}\t{count}\n'
            #frequency_counters_per_sample[sample_barcode][seq] = count
        f.write(counts)
    if not args.extended_output:
        os.remove(files_paths[sample_barcode][FS])

#touch done
with open(done_file, 'w') as f:
    pass

end = time()
logger.info(f'Finished filtering reads. Took {measure_time(start, end)}')


