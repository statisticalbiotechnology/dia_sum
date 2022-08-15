#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:02:02 2022

@author: ptruong
"""

import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['text.usetex'] = True

import numpy as np
import seaborn as sns
import pandas as pd
from parsers.parse_triqler import  parse_triqler
import argparse

def plot_histogram(output, bins = 50):
    fig, axs = plt.subplots(1, 1, figsize=(10,6))
    g = sns.histplot(data=df[df["FDR"] < 1], x="log2FC", hue="specie", kde = False, ax = axs, 
                    bins = np.histogram_bin_edges(df[df["FDR"] < 0.05].log2FC, 
                                                  bins=bins, range=(-2,3), weights=None))
    g.legend_.set_title(None)
    axs.set_xlim([-2, 3])
    axs.set_xlabel("log2 fold change", fontsize=34)
    axs.set_ylabel("Protein count", fontsize=34)
    
    plt.legend(loc='upper right', labels=['Yeast', "HeLa",  r'\textit{E.Coli}'])
    plt.setp(axs.get_legend().get_texts(), fontsize='30')
    
    axs.tick_params(axis='x', which='major', labelsize=32)#labelrotation=90)
    axs.tick_params(axis='y', which='major', labelsize=32)
    
    axs.axvline(x = 0, color = "orange", linestyle = "--", alpha = 0.4)
    axs.axvline(x = -1, color = "green", linestyle = "--", alpha = 0.4)
    axs.axvline(x = 2, color = "blue", linestyle = "--", alpha = 0.4)
    print("Outputting : " + output)
    plt.savefig(output, bbox_inches="tight")

specie_mapper = lambda x: x.split("_")[-1]

def read_triqler(triqler_file):
    df = parse_triqler(triqler_file)
    #df["log2_fold_change"] = -df.log2_fold_change
    df["specie"] = df.protein.map(specie_mapper)
    df.rename({"log2_fold_change":"log2FC", "q_value":"FDR"}, axis = 1, inplace = True)
    df = df[~df.protein.str.contains("DECOY")]
    df.sort_values(by = ["specie"], inplace = True)
    
    df = (df[~((df.iloc[:,df.columns.str.contains("Pedro")] == 1).sum(axis=1) == 6)]) # Remove proteins that are not updates at any sample
    return df

def read_top3(top3_file):
    df = pd.read_csv(top3_file, sep = "\t")
    df.rename({"q":"FDR", "log2(A,B)":"log2FC"}, axis = 1, inplace = True)
    #df.log2FC = -df.log2FC
    df.sort_values(by = ["specie"], inplace = True)
    return df

def read_msstats(msstats_file):
    df = pd.read_csv(msstats_file, sep = ",")
    df.rename({"adj.pvalue":"FDR"}, axis = 1, inplace = True)
    df["specie"] = df.Protein.map(specie_mapper)
    df.sort_values(["specie"], inplace = True)
    df.log2FC = -df.log2FC # added for reverse A/B
    return df

def read_msqrob2(msqrob2_file):
    df = pd.read_csv(msqrob2_file, sep = ",")
    df.rename({"Unnamed: 0":"protein", "logFC":"log2FC", "adjPval":"FDR"}, axis = 1, inplace = True)
    df["specie"] = df.protein.map(specie_mapper)
    df.sort_values(["specie"], inplace = True)
    df.log2FC = -df.log2FC # added for reverse log2FC
    return df



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots differential HeLa vs differential non-HeLa lineplot.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input_file', type=str,
                        help='input file name.')
    
    parser.add_argument('--input_type', type=str,
                        help='the input_type is one of triqler/top3/msstats/msqrob2.')
    
    parser.add_argument('--bins', type=int,
                        help='Number of bins in the histogram.',
                        default=50)
    
    parser.add_argument('--output', type=str,
                        help='output name.')
    
    
    # parse arguments from command line
    args = parser.parse_args()
    input_file = args.input_file
    input_type = args.input_type
    bins = args.bins
    output = args.output
    
    print("Reading in : " + input_file)
    if input_type == "triqler":
        df = read_triqler(input_file)
        plot_histogram(output, bins = bins)
    elif input_type == "top3":
        df = read_top3(input_file)
        plot_histogram(output, bins = bins)
    elif input_type == "msstats":
        df = read_msstats(input_file)
        plot_histogram(output, bins = bins)
    elif input_type == "msqrob2":
        df = read_msqrob2(input_file)
        plot_histogram(output, bins = bins)
    else:
        print("No input_type specified. Exiting!")
    print("Done!")


