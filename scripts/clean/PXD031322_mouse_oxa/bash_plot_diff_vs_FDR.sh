#!/bin/bash

# ST_LT
python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_LT/triqler_results.tsv --top3_input data/ST_LT/top3_results.csv --msstats_input data/ST_LT/msstats_adj_results.csv --msqrob2_input data/ST_LT/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.01 --ylim 400 --output data/ST_LT/diff_protein_vs_FDR_FC_0.415_xlim_0.01.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_LT/triqler_results.tsv --top3_input data/ST_LT/top3_results.csv --msstats_input data/ST_LT/msstats_adj_results.csv --msqrob2_input data/ST_LT/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.05 --ylim 1500 --output data/ST_LT/diff_protein_vs_FDR_FC_0.415_xlim_0.05.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_LT/triqler_results.tsv --top3_input data/ST_LT/top3_results.csv --msstats_input data/ST_LT/msstats_adj_results.csv --msqrob2_input data/ST_LT/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.1 --ylim 1700 --output data/ST_LT/diff_protein_vs_FDR_FC_0.415_xlim_0.1.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_LT/triqler_results.tsv --top3_input data/ST_LT/top3_results.csv --msstats_input data/ST_LT/msstats_adj_results.csv --msqrob2_input data/ST_LT/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.3 --ylim 2200 --output data/ST_LT/diff_protein_vs_FDR_FC_0.415_xlim_0.3.png

# LT_CTRL
python plot_differential_proteins_vs_FDR.py --triqler_input data/LT_CTRL/triqler_results.tsv --top3_input data/LT_CTRL/top3_results.csv --msstats_input data/LT_CTRL/msstats_adj_results.csv --msqrob2_input data/LT_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.01 --ylim 400 --output data/LT_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.01.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/LT_CTRL/triqler_results.tsv --top3_input data/LT_CTRL/top3_results.csv --msstats_input data/LT_CTRL/msstats_adj_results.csv --msqrob2_input data/LT_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.05 --ylim 1500 --output data/LT_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.05.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/LT_CTRL/triqler_results.tsv --top3_input data/LT_CTRL/top3_results.csv --msstats_input data/LT_CTRL/msstats_adj_results.csv --msqrob2_input data/LT_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.1 --ylim 1700 --output data/LT_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.1.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/LT_CTRL/triqler_results.tsv --top3_input data/LT_CTRL/top3_results.csv --msstats_input data/LT_CTRL/msstats_adj_results.csv --msqrob2_input data/LT_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.3 --ylim 2200 --output data/LT_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.3.png

# ST_CTRL
python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_CTRL/triqler_results.tsv --top3_input data/ST_CTRL/top3_results.csv --msstats_input data/ST_CTRL/msstats_adj_results.csv --msqrob2_input data/ST_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.01 --ylim 400 --output data/ST_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.01.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_CTRL/triqler_results.tsv --top3_input data/ST_CTRL/top3_results.csv --msstats_input data/ST_CTRL/msstats_adj_results.csv --msqrob2_input data/ST_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.05 --ylim 1500 --output data/ST_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.05.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_CTRL/triqler_results.tsv --top3_input data/ST_CTRL/top3_results.csv --msstats_input data/ST_CTRL/msstats_adj_results.csv --msqrob2_input data/ST_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.1 --ylim 1700 --output data/ST_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.1.png

python plot_differential_proteins_vs_FDR.py --triqler_input data/ST_CTRL/triqler_results.tsv --top3_input data/ST_CTRL/top3_results.csv --msstats_input data/ST_CTRL/msstats_adj_results.csv --msqrob2_input data/ST_CTRL/msqrob2_results.csv --fc_threshold 0.415 --fdr_threshold 1.00 --xlim 0.3 --ylim 2200 --output data/ST_CTRL/diff_protein_vs_FDR_FC_0.415_xlim_0.3.png

