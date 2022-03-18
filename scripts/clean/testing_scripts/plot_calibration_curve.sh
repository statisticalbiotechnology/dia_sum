#!/bin/bash

python plot_calibration_curve.py --triqler_input ~/data/results/PS/triqler_results/fc_0.48 --top3_input ~/data/results/PS/top3_output_diann.csv --msstats_input ~/data/results/PS/msstat_output.csv --msqrob2_input ~/data/results/PS/msqrob2_results.tsv --fc_threshold 0.48 --output calibration_curve.png

