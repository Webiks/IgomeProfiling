import logging
logger = logging.getLogger('main')

DEFAULT_PARAMS = False

# CONSTANTS:
FILTERED = 'filtered'
INFO = 'info'
FASTQ = 'fastq'
FS = 'fs'
COUNTS = 'counts'

TOTAL_SEQ_AFTER_FILTERS = 'TOTAL SEQ AFTER FILTERS'
DOMAIN_OF_INTEREST = 'recognized_domains'
LEFTOVERS = 'unrecognized_domains'
ALL = 'all_domains'

TOO_MANY_MISTAKES = 'TOO_MANY_MISTAKES'#'TOO_MANY_MISTAKES'
SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE = 'SEQ_WITH_ERRORS_IN_DOMAIN_BARCODE_DOWNSTREAM_SEQUENCE'
WRONG_SEQ_WITH_N = 'WRONG_SEQ_WITH_N'
TOTAL_SEQ = 'TOTAL_SEQ'
TOTAL_SEQS_OK = 'TOTAL_SEQS_OK'

def load_table_to_dict(file_to_read, delimiter ='\t', backwards = 0):
    '''
    parse a delimited file with two columns to a dictionary
    if backwards == 0 keys are taken from the first column o.w., from the second
    open_operator can also be set to gzip.open if the file is zipped
    '''
    d = {}
    values_are_repetative = ''
    with open(file_to_read) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            item1, item2 = line.rstrip().split(delimiter)
            if item1 in d:
                logger.error(f'Cannot load table to dict! key: {item1} appears twice in file. (line: {line.rstrip()})')
                raise KeyError
            if item2 in d:
                values_are_repetative = item2
            d[item1] = item2
    if backwards:
        if values_are_repetative != '':
            logger.error('Cannot load reversed table to dict! value: {} appears twice in file'.format(values_are_repetative))
            #raise ValueError
        d = dict((d[key],key) for key in d)
    return d

def diff(s1, s2):
    assert len(s1) == len(s2)
    return sum([s1[i]!=s2[i] for i in range(len(s1))])

def load_column_to_list(file_path):
    with open(file_path) as f:
        return [line.rstrip() for line in f]

def measure_time(start, end):
    hours = int(end - start) / 3600
    minutes = (int(end - start) % 3600) / 60
    seconds = int(end - start) % 60
    if hours != 0:
        return '%d:%d:%d hours.' % (hours, minutes, seconds)
    elif minutes != 0:
        return '%d:%d minutes.' % (minutes, seconds)
    else:
        return '%d seconds.' % (seconds)

def get_column_from_file(columns_file, column):
    items = []
    with open(columns_file) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            item = line.split('\t')[column]
            items.append(item)
    return items
