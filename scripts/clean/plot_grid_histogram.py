#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:24:48 2022

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

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
    df["log2_fold_change"] = -df.log2_fold_change
    df["specie"] = df.protein.map(specie_mapper)
    df.rename({"log2_fold_change":"log2FC", "q_value":"FDR"}, axis = 1, inplace = True)
    df = df[~df.protein.str.contains("DECOY")]
    df.sort_values(by = ["specie"], inplace = True)
    
    df = (df[~((df.iloc[:,df.columns.str.contains("Pedro")] == 1).sum(axis=1) == 6)]) # Remove proteins that are not updates at any sample
    return df

def read_top3(top3_file):
    df = pd.read_csv(top3_file, sep = "\t")
    df.rename({"q":"FDR", "log2(A,B)":"log2FC"}, axis = 1, inplace = True)
    df.log2FC = -df.log2FC
    df.sort_values(by = ["specie"], inplace = True)
    return df

def read_msstats(msstats_file):
    df = pd.read_csv(msstats_file, sep = ",")
    df.rename({"adj.pvalue":"FDR"}, axis = 1, inplace = True)
    df["specie"] = df.Protein.map(specie_mapper)
    df.sort_values(["specie"], inplace = True)
    #df.log2FC = -df.log2FC # added for reverse A/B
    return df

def read_msqrob2(msqrob2_file):
    df = pd.read_csv(msqrob2_file, sep = ",")
    df.rename({"Unnamed: 0":"protein", "logFC":"log2FC", "adjPval":"FDR"}, axis = 1, inplace = True)
    df["specie"] = df.protein.map(specie_mapper)
    df.sort_values(["specie"], inplace = True)
    #df.log2FC = -df.log2FC # added for reverse log2FC
    return df




def grid_histogram(output, bins = 50, test = False):
    bins = bins
    test = False
    
    triqler_ID = read_triqler("results/ID/triqler_results.csv")
    top3_ID = read_top3("results/ID/top3_results.csv")
    msstats_ID = read_msstats("results/ID/msstats_results.csv")
    msqrob2_ID = read_msqrob2("results/ID/msqrob2_results.csv")
    triqler_PS = read_triqler("results/PS/triqler_results.csv")
    top3_PS = read_top3("results/PS/top3_results.csv")
    msstats_PS = read_msstats("results/PS/msstats_results.csv")
    msqrob2_PS= read_msqrob2("results/PS/msqrob2_results.csv")

    fig2 = plt.figure(constrained_layout=False, figsize=(20,14))
    spec2 = fig2.add_gridspec(ncols=2, nrows=4, wspace=0.0, hspace=0.0)
    
    f2_ax1 = fig2.add_subplot(spec2[0, 0])
    f2_ax2 = fig2.add_subplot(spec2[0, 1])
    f2_ax3 = fig2.add_subplot(spec2[1, 0])
    f2_ax4 = fig2.add_subplot(spec2[1, 1])
    f2_ax5 = fig2.add_subplot(spec2[2, 0])
    f2_ax6 = fig2.add_subplot(spec2[2, 1])
    f2_ax7 = fig2.add_subplot(spec2[3, 0])
    f2_ax8 = fig2.add_subplot(spec2[3, 1])
        
    titles = ["Triqler,ID", "Triqler,PS",
                          "MSqRob2,ID", "MSqRob2,PS",
                          "MSstats,ID", "MSstats,PS",
                          "Top3,ID", "Top3,PS"]
    dfs = [triqler_ID, triqler_PS, 
           msqrob2_ID, msqrob2_PS,
           msstats_ID, msstats_PS,
           top3_ID, top3_PS]
    axs = [f2_ax1, f2_ax2,
           f2_ax3, f2_ax4,
           f2_ax5, f2_ax6,
           f2_ax7, f2_ax8]
    data = list(zip(titles,dfs,axs))
    
    i=0
    for title, df, ax in data:
        g = sns.histplot(data=df[df["FDR"] < 1], x="log2FC", hue="specie", kde = False, ax = ax, 
                     bins = np.histogram_bin_edges(df[df["FDR"] < 0.05].log2FC, 
                                                   bins=bins, range=(-2,3), weights=None))
        ax.set_ylim([0, 999])
        ax.grid(linestyle = ':')
        
        if i == 0:
            ax.set_title("ID", fontsize=28)
        if i == 1:
            ax.set_title("PS", fontsize=28)
        if i % 2 == 1:
            ax.yaxis.set_ticklabels([])
            ax.set_ylabel("")
            ax.tick_params(left=False)
        if i % 2 == 0:
            method = title.split(",")[0]
            ax.set_ylabel(method+"\nProtein count", fontsize=24)
        if i < 6:
            ax.set_xlabel("")
            ax.xaxis.set_ticklabels([])
        if i >= 6:
            ax.set_xlabel("log2(A/B)", fontsize=34)
    
        if test == True:
            ax.set_title(str(title))
        i+=1
         
        ax.tick_params(axis='x', which='major', labelsize=28)#labelrotation=90)
        ax.tick_params(axis='y', which='major', labelsize=28) 
        g.legend_.set_title(None)
        plt.setp(ax.get_legend().get_texts(), fontsize='20')
        
        ax.axvline(x = 0, color = "orange", linestyle = "--", alpha = 0.4)
        ax.axvline(x = -1, color = "green", linestyle = "--", alpha = 0.4)
        ax.axvline(x = 2, color = "blue", linestyle = "--", alpha = 0.4)
        
    print("Outputting : " + output)
    plt.savefig(output, bbox_inches="tight")



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots a grid plot of histograms.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--bins', type=int,
                        help='Number of bins in the histogram.',
                        default=50)
    
    parser.add_argument('--output', type=str,
                        help='output name.')
    
    
    # parse arguments from command line
    args = parser.parse_args()

    bins = args.bins
    output = args.output
    grid_histogram(output, bins)
    print("Plotting gridplot")
    print("Done!")

