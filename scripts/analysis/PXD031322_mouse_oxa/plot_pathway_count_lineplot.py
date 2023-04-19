#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:17:55 2022

@author: ptruong
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

def grid_pathway_plot(count_table_file, percent_DEG = 0, output = "pathway_countplot.png"):
    full_list = pd.read_csv(count_table_file, sep = "\t", index_col = 0)
    full_list.reset_index(inplace = True)
    full_list.rename({"index":"method"}, axis = 1, inplace = True)
    full_list = full_list[~full_list.method.isin(["triqler-intersection", "reported-intersection"])]
    
    
    dfs = []
    titles = []
    percent_DEG = percent_DEG
    for log2FC_fdr in [0.005, 0.01, 0.05]:
        for pathway_fdr in [0.005, 0.01, 0.05]:        
            q = full_list[(full_list["log2FC_fdr"] == log2FC_fdr) &
             (full_list["pathway_fdr"] == pathway_fdr) &
             (full_list["%DEG_in_pathway"] == percent_DEG)
             ]
            dfs.append(q)
            titles.append(f"log2FC_FDR: {log2FC_fdr}, pathway_FDR: {pathway_fdr}, percent_DEG_in_pathway: {percent_DEG}")
    
    
    fig = plt.figure(constrained_layout=False, figsize=(20,14))
    spec = fig.add_gridspec(ncols=3, nrows=3, wspace=0.05, hspace=0.05)
    
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[0, 2])
    ax4 = fig.add_subplot(spec[1, 0])
    ax5 = fig.add_subplot(spec[1, 1])
    ax6 = fig.add_subplot(spec[1, 2])
    ax7 = fig.add_subplot(spec[2, 0])
    ax8 = fig.add_subplot(spec[2, 1])
    ax9 = fig.add_subplot(spec[2, 2])
        
    axs = [ax1, ax2, ax3, 
           ax4, ax5, ax6,
           ax7, ax8, ax9]
    
    data = list(zip(titles,dfs,axs))
    
    i = 0
    for title, df, ax in data:
        g = sns.lineplot(data = df, x = "abs(log2FC)", y = "diff_pathways_Total", hue = "method", ax = ax)
        ax.grid(linestyle = ':')
        #ax.set_ylabel(i)
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_ylim([0,80])
        ax.set_xlim([0.1,0.6])
        ax.tick_params(left=False)
        if i > 0: #Remove legends
                ax.get_legend().remove()
        if i == 0:
            ax.legend(prop={"size":18})
        if i in [0, 3, 6]: # set y-label log2FC FDR
            ax.set_ylabel(title.split(",")[0].split(": ")[1],fontsize=20)
        if i in [6, 7, 8]:
            ax.set_xlabel(title.split(",")[1].split(": ")[-1],fontsize=20)
        # Remove ticks
        if i < 6:
            ax.xaxis.set_ticklabels([])
            ax.tick_params(bottom=False)
        if i in [1,2,4,5,7,8]:
            ax.yaxis.set_ticklabels([])
    
        i+=1
        ax.tick_params(axis='x', which='major', labelsize=14)#labelrotation=90)
        ax.tick_params(axis='y', which='major', labelsize=14) 
        
        if g.legend_ != None:
            g.legend_.set_title(None)
    
    fig.supxlabel("Fold change threeshold\nPathway FDR threshold",fontsize=20)
    fig.supylabel("   Number of pathways\nFold change FDR threshold",fontsize=20)
    fig.suptitle(f"%Pathway is DEG threshold: {percent_DEG}",fontsize=20)
    
    print("Outputting : " + output)
    plt.savefig(output, bbox_inches="tight")

count_table_file = "count_table.tsv"
for percent_DEG in [0, 0.1, 0.2, 0.3, 0.4, 0.5]:
    grid_pathway_plot(count_table_file, percent_DEG = percent_DEG, output = f"pathway_countplot_DEG_{percent_DEG}.png")

