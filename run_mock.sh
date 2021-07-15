#!/bin/bash
rm -rf output
mkdir output && mkdir output/analysis && mkdir output/logs
source .venv/bin/activate
python3 IgOmeProfiling_pipeline.py \
    mock_data/exp12_10M_rows.fastq.gz \
    mock_data/barcode2samplename.txt \
    mock_data/samplename2biologicalcondition.txt \
    output/analysis \
    output/logs
