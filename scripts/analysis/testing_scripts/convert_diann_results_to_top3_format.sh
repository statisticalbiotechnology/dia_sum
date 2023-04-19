#!/bin/bash

python ../top3_formatter/diann.py --input_file /hdd_14T/data/PXD002952/20210614_dataset/result_files_20220214/PS/report_for_top3.tsv --fdr_threshold 0.01 --q_value_column Q.Value --output diann_top3_formatted.csv
