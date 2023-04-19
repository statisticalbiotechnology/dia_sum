#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 16:46:57 2022

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

def get_stats_adjusted(df_inp, step = 0.1):
    df_inp["specie"] = df_inp.protein.map(lambda x:x.split("_")[-1])

    def specie_adj(x):
        if x == "ECOLI":
            return -2
        if x == "YEAST":
            return 1
        else:
            return 0

    df_inp["specie_adj_factor"] = df_inp.specie.map(specie_adj)        
    df_inp["abs(actualFC-estFC)"] = abs(df_inp.log2_fold_change - df_inp["specie_adj_factor"])
    df_inp = df_inp[~df_inp.protein.str.contains("DECOY")]

    means = []
    medians = []
    min_fc_list = []
    max_fc_list = []
    label_name = []
    step = step
    
    for i in np.arange(0,2, step):
        min_fc = i
        max_fc = i+step
        proteins = df_inp[(df_inp["abs(actualFC-estFC)"]>min_fc)&(df_inp["abs(actualFC-estFC)"]<max_fc)].protein
        min_fc_list.append(round(min_fc, 2))
        max_fc_list.append(round(max_fc, 2))
        means.append(df_inp[df_inp.protein.isin(proteins)].num_peptides.mean())
        medians.append(df_inp[df_inp.protein.isin(proteins)].num_peptides.median())
        label_name.append(str(round(min_fc, 2))+"-"+str(round(max_fc,2)))
    res = pd.DataFrame([min_fc_list, max_fc_list, means, medians, label_name], index = ["min_fc", "max_fc", "mean", "median", "label_name"]).T
    res["log2FC"] = (res.min_fc + res.max_fc)/2
    
    
    return res


def plot_number_of_peptides_per_log2FC_range(result_file, output, step = 0.2):
    df_inp = parse_triqler(result_file)
    f, ax = plt.subplots(1, 1, figsize = (24,14))

    res = get_stats_adjusted(df_inp, step = step)
    sns.barplot(data = res, x = "log2FC", y = "mean", 
                color = "tab:blue", ax = ax)
    
    
    ax.legend()
    #ax.axvline(2, linestyle = "--", color="tab:blue", alpha = 0.5)
    #ax.axvline(0, linestyle = "--", color="tab:orange", alpha = 0.5)
    #ax.axvline(-1, linestyle = "--", color="tab:green", alpha = 0.5)
    ax.set_ylim(bottom = 0)
    ax.legend(fontsize=30)

    ax.set_xlabel("Abs(Actual FC - Estimated FC)", fontsize=52)
    ax.set_ylabel("Mean Number of Peptides", fontsize=52)
    
    ax.tick_params(axis='x', which='major', labelsize=36)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=36)
    #ax.legend([r"\textit{E.coli}", "Yeast", "HeLa"], prop = {"size":44})
    

    ax.set_xticklabels(list(res.label_name))
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.2)
    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script takes scatter plot formatted input file and results file to produces a lineplot with number of peptide per fc range.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('--triqler_result_file', type=str,
                        help='Triqler results file.')

    parser.add_argument('--step', type=float, default = 0.2,
                        help='Step to bin the log2FCs. Default: 0.2')

    parser.add_argument('--output', type=str, 
                        help='Output name.')


    # parse arguments from command line
    args = parser.parse_args()
    result_file = args.triqler_result_file
    output = args.output
    step = args.step
    plot_number_of_peptides_per_log2FC_range(result_file, output = output, step = step)
        



    



