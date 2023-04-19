#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 12:12:21 2022

@author: ptruong
"""

import pandas as pd
import numpy as np


def parse_triqler(triqler_output_file):
    """
    Parses triqler output format to pandas dataframe.
    """
    f = open(triqler_output_file, "r")
    lines = f.readlines()
    line = lines.pop(0)
    cols = line.split("\n")[0].split("\t")[:]
    n_cols = len(cols)
    
    data_array = []
    for line in lines:
        line = line.split("\n")[0].split("\t")
        vals = line[:n_cols-1]
        peptides = ";".join(line[n_cols-1:])
        data = vals + [peptides]
        data_array.append(data)
    df = pd.DataFrame(data_array, columns = cols)

    df = pd.concat([df[["protein", "peptides"]], df.drop(["protein", "peptides"], axis = 1).astype(float)], axis = 1)
    
    return df


# 1 - ctrl
# 2 - LT
# 3 - ST
ctrl_lt = parse_triqler("proteins.1vs2.tsv")#.set_index("protein")
ctrl_st = parse_triqler("proteins.1vs3.tsv")#.set_index("protein")
lt_st = parse_triqler("proteins.2vs3.tsv")#.set_index("protein")

ctrl_lt = ctrl_lt[cols]#.rename({"q_value":"ctrl-lt:q_value", "upregulated":"ctrl-lt:upregulated", "log2_fold_change":"ctrl-lt:log2_fold_change"}, axis = 1)
ctrl_st = ctrl_st[cols]#.rename({"q_value":"ctrl-st:q_value", "upregulated":"ctrl-st:upregulated", "log2_fold_change":"ctrl-st:log2_fold_change"}, axis = 1)
lt_st = lt_st[cols]#.rename({"q_value":"lt-st:q_value", "upregulated":"lt-st:upregulated", "log2_fold_change":"lt-st:log2_fold_change"}, axis = 1)

# Make a script to find the clusters 

ctrl_lt["upregulated"] = ctrl_lt.log2_fold_change > 0
ctrl_st["upregulated"] = ctrl_st.log2_fold_change > 0 
lt_st["upregulated"] = lt_st.log2_fold_change > 0

cols = ['protein', 'peptides', 'q_value', 'posterior_error_prob',
       'num_peptides', 'protein_id_posterior_error_prob', 'log2_fold_change', "upregulated"]
ctrl_lt = ctrl_lt[cols]#.rename({"q_value":"ctrl-lt:q_value", "upregulated":"ctrl-lt:upregulated", "log2_fold_change":"ctrl-lt:log2_fold_change"}, axis = 1)
ctrl_st = ctrl_st[cols]#.rename({"q_value":"ctrl-st:q_value", "upregulated":"ctrl-st:upregulated", "log2_fold_change":"ctrl-st:log2_fold_change"}, axis = 1)
lt_st = lt_st[cols]#.rename({"q_value":"lt-st:q_value", "upregulated":"lt-st:upregulated", "log2_fold_change":"lt-st:log2_fold_change"}, axis = 1)

ctrl_lt["condition"] = "ctrl:lt"
ctrl_st["condition"] = "ctrl:st"
lt_st["condition"] = "lt:st"
df = pd.concat([ctrl_lt, ctrl_st, lt_st]).reset_index().drop("index", axis = 1)

#log2fc = df[["ctrl-lt:log2_fold_change", "ctrl-st:log2_fold_change", "lt-st:log2_fold_change"]]

import scipy.stats as stats

stats.f_oneway(df["log2_fold_change"][df["condition"] == "ctrl:lt"],
               df["log2_fold_change"][df["condition"] == "ctrl:st"],
               df["log2_fold_change"][df["condition"] == "lt:st"])


#### TEST

df = pd.read_csv("https://raw.githubusercontent.com/researchpy/Data-sets/master/difficile.csv")
df.drop('person', axis= 1, inplace= True)

# Recoding value from numeric to string
df['dose'].replace({1: 'placebo', 2: 'low', 3: 'high'}, inplace= True)

df.info()




import statsm odels.stats.multicomp as mc

comp = mc.MultiComparison(df["log2_fold_change"], df['condition'])
post_hoc_res = comp.tukeyhsd()
print(post_hoc_res.summary())

post_hoc_res.plot_simultaneous(ylabel= "Condition", xlabel= "Score Difference")



# Manual clustering
cols = ['q_value', 'log2_fold_change', "upregulated"]

ctrl_lt = parse_triqler("proteins.1vs2.tsv").set_index("protein")
ctrl_st = parse_triqler("proteins.1vs3.tsv").set_index("protein")
lt_st = parse_triqler("proteins.2vs3.tsv").set_index("protein")

ctrl_lt["upregulated"] = ctrl_lt.log2_fold_change > 0
ctrl_st["upregulated"] = ctrl_st.log2_fold_change > 0 
lt_st["upregulated"] = lt_st.log2_fold_change > 0


ctrl_lt = ctrl_lt[cols].rename({"q_value":"ctrl-lt:q_value", "upregulated":"ctrl-lt:upregulated", "log2_fold_change":"ctrl-lt:log2_fold_change"}, axis = 1)
ctrl_st = ctrl_st[cols].rename({"q_value":"ctrl-st:q_value", "upregulated":"ctrl-st:upregulated", "log2_fold_change":"ctrl-st:log2_fold_change"}, axis = 1)
lt_st = lt_st[cols].rename({"q_value":"lt-st:q_value", "upregulated":"lt-st:upregulated", "log2_fold_change":"lt-st:log2_fold_change"}, axis = 1)

df = pd.concat([ctrl_lt, ctrl_st, lt_st], axis = 1)

# Clusters 
C1 = ((df["ctrl-lt:upregulated"] == True) & (df["ctrl-st:upregulated"] == True) & (df["lt-st:upregulated"] == False))
C1.sum()
C2 = ((df["ctrl-lt:upregulated"] == False) & (df["ctrl-st:upregulated"] == False) & (df["lt-st:upregulated"] == True))
C2.sum() 
C3 = (df["lt-st:upregulated"] == False)      
C3.sum()
C4 = (df["lt-st:upregulated"] == True)
C4.sum()
C5 = (df["ctrl-st:upregulated"] == True)
C5.sum()
C6 = (df["ctrl-st:upregulated"] == False)
C6.sum()






