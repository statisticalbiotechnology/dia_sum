#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 16:34:47 2022

@author: ptruong
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

def grid_protein_plot(count_table_file, output = "pathway_countplot.png"):
    full_list = pd.read_csv(count_table_file, sep = "\t", index_col = 0)
    full_list.reset_index(inplace = True)
    full_list.rename({"index":"method"}, axis = 1, inplace = True)
    full_list = full_list[~full_list.method.isin(["triqler-intersection", "reported-intersection"])]
    
    
    dfs = []
    titles = []
    percent_DEG = 0
    for log2FC_fdr in [0.005, 0.01, 0.05]:      
            q = full_list[(full_list["log2FC_fdr"] == log2FC_fdr)]
            dfs.append(q)
            titles.append(f"log2FC_FDR: {log2FC_fdr}")
    
    
    fig = plt.figure(constrained_layout=False, figsize=(20,14))
    spec = fig.add_gridspec(ncols=3, nrows=1, wspace=0.05, hspace=0.05)
    
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[0, 2])

        
    axs = [ax1, ax2, ax3]
    
    data = list(zip(titles,dfs,axs))
    
    i = 0
    for title, df, ax in data:
        g = sns.lineplot(data = df, x = "abs(log2FC)", y = "diff_proteins_Total", hue = "method", ax = ax)
        ax.grid(linestyle = ':')
        #ax.set_ylabel(i)
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_ylim([0,2500])
        ax.set_xlim([0.1,0.6])
        ax.tick_params(left=False)
        if i > 0: #Remove legends
                ax.get_legend().remove()
        if i == 0:
            ax.legend(prop={"size":18})
        #if i in [0, 3, 6]: # set y-label log2FC FDR
        #    ax.set_ylabel(title.split(",")[0].split(": ")[1],fontsize=20)
        #if i in [6, 7, 8]:
        ax.set_xlabel(title.split(": ")[-1],fontsize=20)
        # Remove ticks
        #if i < 6:
        #    ax.xaxis.set_ticklabels([])
        #    ax.tick_params(bottom=False)
        if i in [1,2]:
            ax.yaxis.set_ticklabels([])
    
        i+=1
        ax.tick_params(axis='x', which='major', labelsize=14)#labelrotation=90)
        ax.tick_params(axis='y', which='major', labelsize=14) 
        
        if g.legend_ != None:
            g.legend_.set_title(None)
    
    fig.supxlabel("Fold change FDR threshold",fontsize=28)
    fig.supylabel("Number of differential proteins",fontsize=28)
    #fig.suptitle(f"%Pathway is DEG threshold: {percent_DEG}",fontsize=28)
    
    print("Outputting : " + output)
    plt.savefig(output, bbox_inches="tight")

count_table_file = "count_table.tsv"
grid_protein_plot(count_table_file, output = f"protein_countplot.png")
