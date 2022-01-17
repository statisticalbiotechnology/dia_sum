#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 11:03:50 2021

@author: ptruong
"""

import numpy as np
import pandas as pd
import os 
import matplotlib.pyplot as plt

os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811")


df = pd.read_csv("report.tsv", sep = "\t")

# map float in Protein.Ids to check what is the problem and remove these
float_mapper = lambda x : isinstance(x, float)
df["isfloat"] = df["Protein.Ids"].map(float_mapper)
df[df["isfloat"] == True]["Protein.Ids"]
df = df[df["isfloat"] != True]

# get decoy column
decoy_mapper = lambda x:x.split("_")[0] == "DECOY"
df["decoy"] = df["Protein.Ids"].map(decoy_mapper)


#df_run.columns

runs = df.Run.unique()
run = runs[0]


df_run = df[df.Run == run]

df_incorrect = df_run[df_run.decoy == True].sort_values(by = "CScore")
df_correct = df_run[df_run.decoy == False].sort_values(by = "CScore")

df_incorrect["incorrect_q_area"] = df_incorrect["Q.Value"].cumsum()
df_correct["correct_q_area"] =  df_correct["Q.Value"].cumsum()

df_incorrect

##########
df_qVals = pd.DataFrame()

qval_col = "Q.Value"
dfs_qVals = []
for run in runs:
    df_run = df[df.Run == run].sort_values(by = "CScore")
    df_run.set_index("CScore", inplace = True)
    #df_run_correct = df_run[df_run.decoy == False]
    dfs_qVals.append(df_run[qval_col])#.reset_index().drop("index", axis = 1)[qval_col]

for i in dfs_qVals:
    i.plot()
    
    
split_name = lambda x: x.split("_")[5]
df_qVals = df_qVals.rename(split_name, axis = 1)

#(1-df_qVals).plot()
#ax.set_xlabel("ordered PSM by m_Score")
ax.set_xlabel("ordered PSM by CScore")
ax.set_ylabel("m_score")


####### Integrate to FDR #####
# Do we compute FDR on both decoy and target?
# If we look at m_score it should be on the list with decoy and target. We rank all...

# make function from this


def compute_fdr(df_run):
    df_decoy = df_run[df_run.decoy == True]
    df_target = df_run[df_run.decoy == False]
    cscores = []
    fdrs = []
    for i in df_run.index.unique():
        n_target = len(df_target[df_target.index > i])
        n_decoy = len(df_decoy[df_decoy.index > i])
        if (n_target + n_decoy) > 0:
            fdr_i = n_decoy/(n_target + n_decoy)
        fdrs.append(fdr_i)
        cscores.append(i)
    fdr_map = dict(zip(cscores, fdrs))
    df_run["fdr"] = df_run.index.map(fdr_map).fillna(0)
    return df_run

from time import time
start = time()
df_runs = []
for run in runs:
    df_run = df[df.Run == run].sort_values(by = "CScore")
    df_run.set_index("CScore", inplace = True)
    df_run = compute_fdr(df_run)
    df_runs.append(df_run)
    print(time()-start)
end = time()-start
print("Done " + str(end))


eps = 1e-9
df_res = pd.concat(df_runs).reset_index()
df_res.fdr = df_res.fdr + eps #to remove absolute zero and division by zero in searchScore conversion for triqler
df_res.to_csv("report_recomputed_fdr.tsv", sep = "\t", index=False)
     
    
    
    
    
    
    
    
    
    
    