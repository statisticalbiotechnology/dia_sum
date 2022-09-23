#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 16:26:22 2022

@author: ptruong
"""

import pandas as pd

final = pd.read_csv("protein_table_output.tsv", sep = "\t", index_col = 0)



### Scatterplot 

import seaborn as sns
import matplotlib.pyplot as plt

def final_df_to_seaborn_format(final):
    methods = ["triqler", "reported"]
    comparisons = ["CTRL_ST", "CTRL_LT", "LT_ST"]
    
    dfs = []
    #for method in methods:
    for comparison in comparisons:
        triqler_log2fc_col = f"triqler_log2FC_{comparison}"
        triqler_fdr_col = f"triqler_FDR_{comparison}"
        
        reported_log2fc_col = f"reported_log2FC_{comparison}"
        reported_fdr_col = f"reported_FDR_{comparison}"

        df = final[[triqler_log2fc_col, triqler_fdr_col,
                    reported_log2fc_col, reported_fdr_col]].rename({f"triqler_log2FC_{comparison}": "triqler_log2FC", 
                                                  f"triqler_FDR_{comparison}": "triqler_FDR",
                                                  f"reported_log2FC_{comparison}": "reported_log2FC",
                                                  f"reported_FDR_{comparison}": "reported_FDR"}, axis = 1)
        #df["method"] = method
        df["comparison"] = comparison
        df.reset_index(inplace = True)
        dfs.append(df)

    res = pd.concat(dfs)
    return res

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes



df = final_df_to_seaborn_format(final)

f, ax = plt.subplots(1, 1, figsize = (15,10))

sns.scatterplot(data=df, x="triqler_log2FC", y="reported_log2FC", hue = "comparison", ax = ax)

add_identity(ax, color='k', ls='--')

