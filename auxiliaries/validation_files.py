from auxiliaries.pipeline_auxiliaries import load_table_to_dict
import re
import datetime


def is_same_samples(samples2bc_dict, barcode2samples_dict, samplename2biologicalcondition_path, barcode2samplename_path, logger):
    # verify that file sample2biologicalcondition and file barcode2sample have the same samples.
    samples_from_sample2bc = sorted(samples2bc_dict.keys())
    samples_from_barcode2sample = sorted(barcode2samples_dict.values())
    if not samples_from_sample2bc == samples_from_barcode2sample:
        barcodes = set(samples_from_barcode2sample)
        bcs = set(samples_from_sample2bc)
        diff_barcodes = barcodes - bcs
        diff_bcs = bcs - barcodes
        logger.error(f'The files: {samplename2biologicalcondition_path} and {barcode2samplename_path} not contain the same samples!!')
        if len(diff_barcodes):
            logger.error(f"Barcodes has samples which doesn't appear in biological conditions: {diff_barcodes}")
        if len(diff_bcs):
            logger.error(f"Biological conditions has samples which doesn't appear in barcodes: {diff_bcs}")
        return False
    return True


def is_valid_data(path_file ,dict_data, logger):
    # verify the structure of the file -  should look like: word\tword\n
    if  not dict_data:
        logger.info(f'There is not data in file: {path_file}')
        return False
    regex = "^[A-Za-z0-9_]+$"
    for key, value in  dict_data.items():
        match_key = re.match(regex, key)
        match_value = re.match(regex, value)
        if match_key is None or match_value is None:
            logger.info(f'The structure of file: {path_file} is not valid')
            return False
    return True        


def load_and_validate_input_table(path, error, logger):
    is_valid = True
    data = None
    try:
        data =  load_table_to_dict(path, error, is_validate_json=True)
    except:
        logger.error(f'{datetime.datetime.now()}: Can not load the file - {path}')
        is_valid = False
    if is_valid:
        is_valid = is_valid_data(path, data, logger) 
    return is_valid, data


def is_input_files_valid(samplename2biologicalcondition_path, barcode2samplename_path, logger):
    samples2bc_valid = True
    barcode2samples_valid = True
    barcode2sample_data = {}
    sample2bc_data = {}

    if barcode2samplename_path:
        barcode2samples_valid, barcode2sample_data = load_and_validate_input_table(barcode2samplename_path, \
            'Barcode {} belongs to more than one sample!!', logger)
    
    if samplename2biologicalcondition_path:
        samples2bc_valid, sample2bc_data = load_and_validate_input_table(samplename2biologicalcondition_path, \
            'Samples {} belongs to more than one bc!!', logger)
    
    if not samples2bc_valid or not barcode2samples_valid:
        return False
    if barcode2samplename_path and samplename2biologicalcondition_path:
        return is_same_samples(sample2bc_data, barcode2sample_data, samplename2biologicalcondition_path, barcode2samplename_path, logger)
    return True     


if __name__ == '__main__':
    import sys
    from logging import basicConfig, getLogger, INFO
    basicConfig(level=INFO)
    logger = getLogger('main')

    barcode2samples_path = sys.argv[1] if len(sys.argv) > 1 else ""
    samples2bc_path = sys.argv[2] if len(sys.argv) > 2 else ""
    is_valid = is_input_files_valid(samples2bc_path, barcode2samples_path, logger)
    print(f'bar2sam="{barcode2samples_path}", sam2bc="{samples2bc_path}", is valid: {is_valid}')
