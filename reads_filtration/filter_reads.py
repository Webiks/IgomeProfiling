import datetime
import gzip
import os
import sys
if os.path.exists('/groups/pupko/orenavr2/'):
    src_dir = '/groups/pupko/orenavr2/igomeProfilingPipeline/src'
elif os.path.exists('/Users/Oren/Dropbox/Projects/'):
    src_dir = '/Users/Oren/Dropbox/Projects/gershoni/src'
else:
    src_dir = '.'
sys.path.insert(0, src_dir)

# needs $src_dir in path
from auxiliaries.pipeline_auxiliaries import load_table_to_dict, fail


def get_barcodes_dictionaries(barcode_to_samplename, output_dir, gz) -> {str: {str: str}}:
    barcode2filepaths = {}
    barcode2filehandlers = {}
    barcode2statistics = {}
    barcode2info_filenames = {}

    open_function = open
    mode = 'w'
    extension = ''
    if gz:
        open_function = gzip.open
        mode = 'wt'
        extension = '.gz'

    for barcode in barcode_to_samplename:
        sample_name = barcode_to_samplename[barcode]  # f'{barcode}_{os.path.split(fastq_file)[-1]}'

        # create a sub dir for that barcode (under output_dir)
        subdir_path = f'{output_dir}/{sample_name}'
        os.makedirs(subdir_path, exist_ok=True)

        # add filenames
        barcode2filepaths[barcode] = {}
        # translated seq
        barcode2filepaths[barcode]['fastq'] = f'{subdir_path}/{sample_name}.fastq{extension}'
        # translated seq
        barcode2filepaths[barcode]['fna'] = f'{subdir_path}/{sample_name}.fna{extension}'
        # translated seq
        barcode2filepaths[barcode]['faa'] = f'{subdir_path}/{sample_name}.faa{extension}'
        # where in the filter seqs dropped
        barcode2filepaths[barcode]['filtration_log'] = f'{subdir_path}/{sample_name}.filtration_log.txt{extension}'

        # add filehandlers
        barcode2filehandlers[barcode] = {}
        barcode2filehandlers[barcode]['fastq'] = open_function(barcode2filepaths[barcode]['fastq'], mode)
        barcode2filehandlers[barcode]['fna'] = open_function(barcode2filepaths[barcode]['fna'], mode)
        barcode2filehandlers[barcode]['faa'] = open_function(barcode2filepaths[barcode]['faa'], mode)
        barcode2filehandlers[barcode]['filtration_log'] = open_function(barcode2filepaths[barcode]['filtration_log'], mode)

        barcode2statistics[barcode] = dict.fromkeys(['legal_barcode',
                                     'poor_quality_barcode',
                                     'high_quality_barcode',
                                     'too_many_mistakes',
                                     'stop_codon',
                                     'too_short',
                                     'total_translated_sequences',
                                     'uag'], 0)
        barcode2statistics[barcode]['lib_type'] = {}

        barcode2info_filenames[barcode] = f'{subdir_path}/{sample_name}_info.txt'

    return barcode2filepaths, barcode2filehandlers, barcode2statistics, barcode2info_filenames


def str_diff(s1, s2):
    # s1 should be equal or shorter. If not, swap.
    if len(s2) < len(s1):
        s1, s2 = s2, s1
    diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    return diff


def write_header(f_handler, txt):
    f_handler.write('_' * 80 + f'\n{txt}')


def filter_reads(argv, fastq_path, parsed_fastq_results, logs_dir,
                 done_path, barcode2samplename_path, name_summary_file,
                 left_construct, right_construct, max_mismatches_allowed,
                 min_sequencing_quality, minimal_length_required, gz):

    start_time = datetime.datetime.now()
    from auxiliaries.pipeline_auxiliaries import nnk_table

    left_construct_length = len(left_construct)
    right_construct_length = len(right_construct)

    logger.info(f'{datetime.datetime.now()}: Loading barcode2samplename file from:\n{barcode2samplename_path}')

    if os.path.isdir(fastq_path):
        fastq_files = [os.path.join(fastq_path, file) for file in os.listdir(fastq_path) if file.endswith('.fastq.gz')]
        total_files = len(fastq_files)
        logger.info(f'{total_files} fastq files found: {fastq_files}')
    else:
        fastq_files = [fastq_path]
        total_files = 1

    barcode2samplename = load_table_to_dict(barcode2samplename_path, 'Barcode {} belongs to more than one sample!!')
    assert len(barcode2samplename) > 0, f'No barcodes were found in {barcode2samplename_path}'  # TODO: add informative error to log

    # get barcode's length by taking the length of an arbitrary barcode
    for barcode in barcode2samplename:
        barcode_len = len(barcode)
        break

    # initialize dictionaries that map barcode to output file names, file handlers and statistics
    barcode2filepaths, barcode2filehandlers, barcode2statistics, barcode2info_filenames = get_barcodes_dictionaries(barcode2samplename, parsed_fastq_results, gz)

    lib_types = set()  # keeps all the library types seen
    record_num = 0  # so the line number will represent the actual sequence inside each record
    for i, fastq_file in enumerate(fastq_files):
        logger.info(f'Reading and filtering {i+1}/{total_files}: {fastq_file}')
        with gzip.open(fastq_file, 'rb') as fastq_f:
            for line1 in fastq_f:
                # ignore first line of each record
                line1 = line1.decode("utf-8").rstrip()
                # second line contains the read itself
                dna_read = fastq_f.readline().decode("utf-8").rstrip()
                # ignore third line of each record
                line3 = fastq_f.readline().decode("utf-8").rstrip()
                # fourth line contains the sequencing quality
                quality = fastq_f.readline().decode("utf-8").rstrip()


                if record_num % 1_000_000 == 0:
                    logger.info(f'{datetime.datetime.now()}: {record_num} records have been processed...')

                if not line1.startswith('@'):
                    raise TypeError('NGS file has an invalid format')

                record_num += 1

                # check that barcode exists
                barcode = dna_read[:barcode_len]
                if barcode not in barcode2samplename:
                    continue

                barcode2statistics[barcode]['legal_barcode'] += 1

                # check that the barcode quality is above the required threshold
                barcode_quality = quality[:barcode_len]
                lowest_scored_character = min([ord(c) for c in barcode_quality])
                if lowest_scored_character < min_sequencing_quality:
                    barcode2statistics[barcode]['poor_quality_barcode'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"{barcode2statistics['legal_barcode']}\tlow quality barcode ({lowest_scored_character})\t{dna_read}\n")
                    continue
                barcode2statistics[barcode]['high_quality_barcode'] += 1

                # verify left construct
                current_left_construct = dna_read[barcode_len: barcode_len + left_construct_length]
                # if regex.match(f'({current_left_construct})' + '{s<='f"{max_mismatches_allowed}"'}', left_construct) == None:
                num_of_mismatches_on_the_left = str_diff(current_left_construct, left_construct)
                if num_of_mismatches_on_the_left > max_mismatches_allowed:
                    barcode2statistics[barcode]['too_many_mistakes'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"Sequence number {barcode2statistics[barcode]['legal_barcode']}\tleft construct has {num_of_mismatches_on_the_left} (max mismatches allowed is: {max_mismatches_allowed})\t{current_left_construct}\n")
                    continue

                # extract random fragment (according to the currently assumed library)
                random_dna_start = barcode_len + left_construct_length
                rest_of_read = dna_read[random_dna_start:]

                random_peptide = ''
                has_stop_codon = False
                q_in_peptide = False
                i = 0
                for i in range(0, len(rest_of_read), 3):
                    codon = rest_of_read[i: i+3]
                    if codon == "TGA" or codon == "TAA":
                        has_stop_codon = True
                        barcode2statistics[barcode]['stop_codon'] += 1
                        barcode2filehandlers[barcode]['filtration_log'].write(f"Sequence number {barcode2statistics[barcode]['legal_barcode']}\tcontains a stop codon\t{rest_of_read[i:i+3]}\n")
                        break
                    if len(codon) < 3:
                        break
                    if codon[0] not in 'ACGT' or codon[1] not in 'ACGT':  # unrecognized dna base, e.g., N
                        break
                    if codon[2] not in 'GT':  # K group (of NNK) = G or T
                        break
                    random_peptide += nnk_table[codon]
                    if codon == 'TAG':
                        q_in_peptide = True

                if has_stop_codon:
                    # all set and documented. We can continue to the next read...
                    continue

                random_dna_end = i  # where we stopped seeing NNK
                rest_of_read = rest_of_read[random_dna_end:] # current right construct + dna remnants

                # if the right construct is too short, just try to match what it has...
                num_of_mismatches_on_the_right = str_diff(rest_of_read, right_construct)
                if num_of_mismatches_on_the_right > max_mismatches_allowed - num_of_mismatches_on_the_left:
                    barcode2statistics[barcode]['too_many_mistakes'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"Sequence number {barcode2statistics[barcode]['legal_barcode']}\tflanking constructs have {num_of_mismatches_on_the_left + num_of_mismatches_on_the_right} (max mismatches allowed is: {max_mismatches_allowed})\t{current_left_construct} {rest_of_read[:right_construct_length]}\n")
                    # all set and documented. We can continue to the next read...
                    continue

                if (len(random_peptide) < minimal_length_required or
                    (random_peptide.startswith('C') and random_peptide.endswith('C') and
                    len(random_peptide)-2 < minimal_length_required)):  # minimum required length (excluding flanking Cysteine)
                    barcode2statistics[barcode]['too_short'] += 1
                    barcode2filehandlers[barcode]['filtration_log'].write(f"Sequence number {barcode2statistics[barcode]['legal_barcode']}\trandom peptide is too short\t{random_peptide}\n")
                    # all set and documented. We can continue to the next read...
                    continue

                # reached here? read is valid!! woohoo
                barcode2statistics[barcode]['total_translated_sequences'] += 1
                if q_in_peptide:
                    barcode2statistics[barcode]['uag'] += 1

                # update library counts
                current_lib_type = f'{len(random_peptide)}'
                if random_peptide.startswith('C') and random_peptide.endswith('C'):
                    # A random peptide from a cys loop library
                    current_lib_type = f'C{len(random_peptide)-2}C'
                if current_lib_type not in lib_types:
                    lib_types.add(current_lib_type)
                barcode2statistics[barcode]['lib_type'][current_lib_type] = barcode2statistics[barcode]['lib_type'].get(current_lib_type, 0) + 1

                # add record to a dedicated fastq file
                barcode2filehandlers[barcode]['fastq'].write(f'{line1}\n{dna_read}\n{line3}\n{quality}\n')
                # add dna read to a dedicated fna file
                barcode2filehandlers[barcode]['fna'].write(
                    f">Seq_{barcode2statistics[barcode]['legal_barcode']}_Lib_{current_lib_type}\n{dna_read[random_dna_start: random_dna_end]}\n")
                # add peptide (translated dna) to a dedicated faa file
                barcode2filehandlers[barcode]['faa'].write(
                    f">Seq_{barcode2statistics[barcode]['legal_barcode']}_Lib_{current_lib_type}\n{random_peptide}\n")

    for barcode in barcode2samplename:
        barcode2filehandlers[barcode]['fastq'].close()
        barcode2filehandlers[barcode]['faa'].close()
        barcode2filehandlers[barcode]['fna'].close()
        barcode2filehandlers[barcode]['filtration_log'].close()

    lib_types = sorted(lib_types, key=lambda x: int(x) if x.isdigit() else 100 + int(x[1:-1]))  # no C's then C's

    for barcode in barcode2samplename:
        with open(barcode2info_filenames[barcode], 'w') as f:
            f.write(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}\n')
            f.write(f'reads_filtration function is invoked with the following parameters:\n')
            f.write(f'fastq_path = {fastq_path}\n'
                    f'parsed_fastq_results = {parsed_fastq_results}\n'
                    f'logs_dir = {logs_dir}\n'
                    f'barcode2samplename_path = {barcode2samplename_path}\n'
                    f'left_construct = {left_construct}\n'
                    f'right_construct = {right_construct}\n'
                    f'max_mismatches_allowed = {max_mismatches_allowed}\n'
                    f'min_sequencing_quality = {min_sequencing_quality}\n')

            total_filtered_reads_per_barcode = barcode2statistics[barcode]['poor_quality_barcode'] + \
                                   barcode2statistics[barcode]['too_many_mistakes'] + \
                                   barcode2statistics[barcode]['stop_codon'] + \
                                   barcode2statistics[barcode]['too_short']

            sanity_check = barcode2statistics[barcode]['poor_quality_barcode'] + barcode2statistics[barcode]['high_quality_barcode'] == barcode2statistics[barcode]['legal_barcode']
            assert sanity_check, 'poor_quality + high_quality != total'

            write_header(f, f"Total sequences with a appropriate barcode: {barcode2statistics[barcode]['legal_barcode']}\n")
            write_header(f, f"Total number of reads that were filtered (due to the following reasons) -> {total_filtered_reads_per_barcode}\n")
            f.write(f"Poor quality barcode -> {barcode2statistics[barcode]['poor_quality_barcode']}\n" 
                    f"More than {max_mismatches_allowed} mistakes in the flanking constant sequences -> {barcode2statistics[barcode]['too_many_mistakes']}\n"
                    f"Nonsense stop codon -> {barcode2statistics[barcode]['stop_codon']}\n"           
                    f"Too short NNK sequences -> {barcode2statistics[barcode]['too_short']}\n")

            write_header(f, f"\n\nTOTAL NUMBER OF TRANSLATED SEQUENCES -> {barcode2statistics[barcode]['total_translated_sequences']}\n\n")

            write_header(f, 'Translated sequences (per library):\n')
            for lib_type in lib_types:
                f.write(f"{lib_type} -> {barcode2statistics[barcode]['lib_type'].get(lib_type, 0)}\n")

            write_header(f, f"Total number of dna sequences with UAG (q) -> {barcode2statistics[barcode]['uag']}\n")


    logger.info('Writing execution summary file...')

    with open(f'{parsed_fastq_results}/{name_summary_file}', 'w') as log_f:
        log_f.write(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}\n')
        log_f.write(f'reads_filtration function is invoked with the following parameters:\n')
        log_f.write(f'fastq_path = {fastq_path}\n'
                    f'parsed_fastq_results = {parsed_fastq_results}\n'
                    f'logs_dir = {logs_dir}\n'
                    f'barcode2samplename_path = {barcode2samplename_path}\n'
                    f'left_construct = {left_construct}\n'
                    f'right_construct = {right_construct}\n'
                    f'max_mismatches_allowed = {max_mismatches_allowed}\n'
                    f'min_sequencing_quality = {min_sequencing_quality}\n')
        log_f.write(f'Barcodes = {" ".join(barcode2samplename)}\n'
                    f'Sample names = {" ".join(barcode2samplename[barcode] for barcode in barcode2samplename)}\n\n')

        reads_with_legal_barcode = sum(barcode2statistics[barcode]['legal_barcode'] for barcode in barcode2samplename)
        write_header(log_f, f'Total number of dna sequences (i.e., fastq records) processed -> {record_num}\n')
        write_header(log_f, f'Total number of dna sequences with an INAPPROPRIATE barcode -> {record_num - reads_with_legal_barcode}\n')
        write_header(log_f, f'Total number of dna sequences with a LEGAL barcode -> {reads_with_legal_barcode}\n')

        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["legal_barcode"]}\n')

        total_filtered_reads = sum(barcode2statistics[barcode]['poor_quality_barcode'] +
                               barcode2statistics[barcode]['too_many_mistakes'] +
                               barcode2statistics[barcode]['stop_codon'] +
                               barcode2statistics[barcode]['too_short'] for barcode in barcode2samplename)

        write_header(log_f, f"Total number of reads that were filtered (due to the following reasons) -> {total_filtered_reads}\n")

        log_f.write(f"Poor quality barcode -> {sum(barcode2statistics[barcode]['poor_quality_barcode'] for barcode in barcode2samplename)}\n"
                    f"More than {max_mismatches_allowed} mistakes in the flanking constant sequences -> "
                               f"{sum(barcode2statistics[barcode]['too_many_mistakes'] for barcode in barcode2samplename)}\n"
                    f"Nonsense stop codon -> {sum(barcode2statistics[barcode]['stop_codon'] for barcode in barcode2samplename)}\n"
                    f"Not NNK sequences -> {sum(barcode2statistics[barcode]['too_short'] for barcode in barcode2samplename)}\n")

        total_translated_sequences = sum(barcode2statistics[barcode]['total_translated_sequences'] for barcode in barcode2samplename)
        write_header(log_f, f'\n\nTOTAL NUMBER OF TRANSLATED SEQUENCES -> {total_translated_sequences}\n\n')

        write_header(log_f, 'Translated sequences (per library):\n')
        for lib_type in lib_types:
            log_f.write(f"{lib_type} -> {sum(barcode2statistics[barcode]['lib_type'].get(lib_type, 0) for barcode in barcode2samplename)}\n")

        write_header(log_f, 'Translated peptides (per barcode):\n')
        for barcode in barcode2samplename:
            log_f.write(f'{barcode} -> {barcode2statistics[barcode]["total_translated_sequences"]}\n')

        uag = sum(barcode2statistics[barcode]["uag"] for barcode in barcode2samplename)
        write_header(log_f, f'Total number of dna sequences with UAG (q) -> {uag}\n')

        end_time = datetime.datetime.now()
        write_header(log_f, f'Total running time: {str(end_time-start_time)[:-3]}')

        log_f.close()

        logger.info(f'Started processing at {start_time}')
        logger.info(f'Done Processing! at {end_time}')
        logger.info(f'Total running time: {str(end_time-start_time)[:-3]}')

    with open(done_path, 'w') as f:
        f.write(" ".join(argv) + '\n')


if __name__ == '__main__':
    print(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq_path', type=str, help='A fastq file to parse')
    parser.add_argument('parsed_fastq_results', type=str, help='folder output')
    parser.add_argument('logs_dir', type=str, help='logs folder')
    parser.add_argument('done_file_path', help='A path to a file that signals that the script finished running successfully.')
    parser.add_argument('barcode2samplename', type=str, help='A path to the barcode to sample name file')
    parser.add_argument('--name_summary_file', default='summary_log.txt', type=str, help='A name for summary file of filter reads.')
    parser.add_argument('--error_path', type=str, help='a file in which errors will be written to')
    parser.add_argument('--left_construct', type=str, default="CAACGTGGC", help='left constant sequence')
    parser.add_argument('--right_construct', type=str, default="GCCT", help='right constant sequence')
    parser.add_argument('--max_mismatches_allowed', type=int, default=1,
                        help='number of mismatches allowed together in both constructs')
    parser.add_argument('--min_sequencing_quality', type=int, default=38,
                        help='Minimum average sequencing threshold allowed after filtration'
                             'for more details, see: https://en.wikipedia.org/wiki/Phred_quality_score')
    # more than 12 aa-long random peptide (e.g., C12C) is irrelevant for 51bps-long reads
    # parser.add_argument('--lib_types', type=str.upper, default='6,C6C,8,C8C,10,C10C,12', help='OBSOLETE: Ignore this param. CxC,x')
    parser.add_argument('--minimal_length_required', default=3, type=int,
                        help='Shorter peptides will be discarded')
    parser.add_argument('--gz', action='store_true', help='gzip fastq, filtration_log, fna, and faa files')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase output verbosity')
    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('main')

    error_path = args.error_path if args.error_path else os.path.join(args.parsed_fastq_results, 'error.txt')

    filter_reads(sys.argv, args.fastq_path, args.parsed_fastq_results, args.logs_dir,
                 args.done_file_path, args.barcode2samplename, args.name_summary_file,
                 args.left_construct, args.right_construct, args.max_mismatches_allowed,
                 args.min_sequencing_quality, args.minimal_length_required, args.gz)
