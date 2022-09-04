#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:36:52 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from parsers.parse_triqler import parse_triqler
import argparse
#sns.set_context("talk")
from matplotlib import rcParams
rcParams['text.usetex'] = True


def get_stats(df, df_inp, step = 0.2):
    """
    use df to select proteins. They should be specie specific since each 
    density is its own specie in the histogram.
    """
    means = []
    medians = []
    min_fc_list = []
    max_fc_list = []
    step = step
    for i in np.arange(-3, 3, step):
        min_fc = i
        max_fc = i+step
        proteins = df[(df["log2(A,B)"] > min_fc) & (df["log2(A,B)"] < max_fc)].ProteinName
        min_fc_list.append(round(min_fc, 2))
        max_fc_list.append(round(max_fc, 2))
        means.append(df_inp[df_inp.protein.isin(proteins)].num_peptides.mean())
        medians.append(df_inp[df_inp.protein.isin(proteins)].num_peptides.median())
    res = pd.DataFrame([min_fc_list, max_fc_list, means, medians], index = ["min_fc", "max_fc", "mean", "median"]).T
    res["log2FC"] = (res.min_fc + res.max_fc)/2
    return res

def plot_number_of_peptides_per_log2FC_range(result_file, input_file, step = 0.2):
    df = pd.read_csv(input_file, sep = "\t")
    df_inp = parse_triqler(result_file)
    step = step
    f, ax = plt.subplots(1, 1, figsize = (24,14))

    res_ecoli = get_stats(df[df.specie == "ECOLI"], df_inp, step = step) # g
    res_yeast = get_stats(df[df.specie == "YEAST"], df_inp, step = step) # o
    res_human = get_stats(df[df.specie == "HUMAN"], df_inp, step = step) # b
    sns.lineplot(data = res_ecoli, x = "log2FC", y = "median", linestyle='-', color = "tab:blue", ax = ax, label = r'\textit{E.coli}')
    sns.lineplot(data = res_yeast, x = "log2FC", y = "median", linestyle='-', color = "tab:green", ax = ax, label = "Yeast")
    sns.lineplot(data = res_human, x = "log2FC", y = "median", linestyle='-', color = "tab:orange", ax = ax, label = "HeLa")
    #plt.legend(labels=["ECOLI","YEAST", "HUMAN"])
    
    
    
    
    ax.legend()
    ax.axvline(2, linestyle = "--", color="tab:blue", alpha = 0.5)
    ax.axvline(0, linestyle = "--", color="tab:orange", alpha = 0.5)
    ax.axvline(-1, linestyle = "--", color="tab:green", alpha = 0.5)
    ax.set_ylim(bottom = 0)
    ax.legend(fontsize=30)

    ax.set_xlabel("Log2(A/B)", fontsize=34)
    ax.set_ylabel("Median Number of Peptides", fontsize=34)
    
    ax.tick_params(axis='x', which='major', labelsize=32)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=32)
    
    
    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script takes scatter plot formatted input file and results file to produces a lineplot with number of peptide per fc range.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--triqler_input_file', type=str,
                        help='Scatterplot converter input file from triqler.')

    parser.add_argument('--triqler_result_file', type=str,
                        help='Triqler results file.')

    parser.add_argument('--step', type=float, default = 0.2,
                        help='Step to bin the log2FCs. Default: 0.2')

    parser.add_argument('--output', type=str, 
                        help='Output name.')


    # parse arguments from command line
    args = parser.parse_args()
    input_file = args.triqler_input_file
    result_file = args.triqler_result_file
    output = args.output
    step = args.step
    df = pd.read_csv(input_file, sep = "\t")
    plot_number_of_peptides_per_log2FC_range(result_file, input_file, step = step)
        



    


