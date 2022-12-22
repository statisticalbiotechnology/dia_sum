#!/bin/bash

python plot_differential_proteins_vs_FDR.py --triqler_input data/triqler_results.tsv --top3_input data/top3_results.csv --msstats_input data/msstats_results.csv --msqrob2_input data/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --output diff_protein_vs_FDR_FC_0.415.png
