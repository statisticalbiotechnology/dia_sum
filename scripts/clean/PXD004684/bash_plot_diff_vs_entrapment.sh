#!/bin/bash

python plot_differential_protein_vs_entrapment.py --triqler_input data/triqler_results.tsv --top3_input data/top3_results.csv --msstats_input data/msstats_results.csv --msqrob2_input data/msqrob2_results.csv --fc_threshold 0.92 --output data/diff_entrapment_fc_0.92.png

python plot_differential_protein_vs_entrapment.py --triqler_input data/triqler_results.tsv --top3_input data/top3_results.csv --msstats_input data/msstats_results.csv --msqrob2_input data/msqrob2_results.csv --fc_threshold 1.00 --output data/diff_entrapment_fc_1.00.png

python plot_differential_protein_vs_entrapment.py --triqler_input data/triqler_results.tsv --top3_input data/top3_results.csv --msstats_input data/msstats_results.csv --msqrob2_input data/msqrob2_results.csv --fc_threshold 1.4 --output data/diff_entrapment_fc_1.4.png

python plot_differential_protein_vs_entrapment.py --triqler_input data/triqler_results.tsv --top3_input data/top3_results.csv --msstats_input data/msstats_results.csv --msqrob2_input data/msqrob2_results.csv --fc_threshold 1.6 --output data/diff_entrapment_fc_1.6.png