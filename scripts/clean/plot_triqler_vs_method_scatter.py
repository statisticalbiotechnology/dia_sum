#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:22:40 2022

@author: ptruong
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from matplotlib import rcParams

sns.set_context("talk")
rcParams["text.usetex"] = True

#df_top3 = pd.read_csv("top3_results.csv", sep = "\t")
#df = pd.read_csv("msstats_scatter_format_input.csv", sep = "\t")
#df = pd.read_csv("msqrob2_scatter_format_input.csv", sep = "\t")
#df = pd.read_csv("triqler_scatter_format_input.csv", sep = "\t")
#df["log2(A,B)"] = -df["log2(A,B)"]
#df["log2(A,B)"] = df["2"]-df["1"]


def plot_scatterplot(df_triqler, df_method, output, fdr_threshold = 1.00):
    
    df_triqler = pd.read_csv(df_triqler, sep = "\t")
    df_method = pd.read_csv(df_method, sep = "\t")
    df_triqler = df_triqler[df_triqler["q"] < fdr_threshold]
    df_method = df_method[df_method["q"] < fdr_threshold]
    df_triqler.set_index("ProteinName", inplace = True)
    df_method.set_index("ProteinName", inplace = True)
    df_triqler["triqler"] = df_triqler["log2(A,B)"]
    df_method["method"] = df_method["log2(A,B)"]
    
    df = pd.concat([df_triqler.triqler, df_method.method], axis = 1).dropna()
    df = df[np.isfinite(df["triqler"])]
    df = df[np.isfinite(df["method"])]
    df.reset_index(inplace = True)
    df["specie"] = df.ProteinName.map(lambda x:x.split("_")[-1])
    
    f, ax = plt.subplots(1, 1, figsize = (15,10))
    sns.regplot(data = df[df.specie == "ECOLI"], x = "triqler", y = "method", ax = ax, line_kws = {"ls":"--"}, label ="ECOLI", fit_reg = False)
    sns.regplot(data = df[df.specie == "HUMAN"], x = "triqler", y = "method", ax = ax, line_kws = {"ls":"--"}, label = "HUMAN", fit_reg = False)
    sns.regplot(data = df[df.specie == "YEAST"], x = "triqler", y = "method", ax = ax, line_kws = {"ls":"--"}, label = "YEAST", fit_reg = False)
    ax.legend()
    ax.axhline(2, linestyle = "--", color="tab:blue", alpha = 0.5)
    ax.axhline(0, linestyle = "--", color="tab:orange", alpha = 0.5)
    ax.axhline(-1, linestyle = "--", color="tab:green", alpha = 0.5)

    ax.set_ylim([-2,3])
    #ax.set_xlim([0,10])
    ax.set_xlabel("Log2(B)", fontsize=38)
    ax.set_ylabel("Log2(A/B)", fontsize=38)
    
    ax.tick_params(axis='x', which='major', labelsize=38)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=38)
    ax.legend(["Yeast", "HeLa", r"\textit{E.coli}"], prop = {"size":40}) 
    
    #plt.title("std/mu ratio for log-transformed peptide values")
        
    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script takes scatter plot formatted input file and produces a scatter plot.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--triqler_input', type=str,
                        help='Scatterplot converter triqler input file.')

    parser.add_argument('--method_input', type=str,
                        help='Scatterplot converter top3, msstats, msqrob2 input file.')
    
    parser.add_argument('--output', type=str,
                        help='Output name.')

    parser.add_argument('--fdr_threshold', type=float, default = 1.00,
                        help='FDR threshold limit. Default: 1.00')


    # parse arguments from command line
    args = parser.parse_args()
    triqler_input = args.triqler_input
    method_input = args.method_input
    fdr_threshold = args.fdr_threshold
    output = args.output
    plot_scatterplot(triqler_input, method_input, output, fdr_threshold)
    
    
