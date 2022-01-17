#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 15:30:31 2022

@author: ptruong
"""

# diann to msqrob2 parser

import os
import pandas as pd
import numpy as np
os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811")

fdr_treshold = 0.01

df = pd.read_csv("report_recomputed_fdr.tsv", sep = "\t")
df.sort_values("CScore", ascending = False, inplace = True)
df = df.groupby("Stripped.Sequence").head(1) # select top scoring PSM

# We need to recompute fdr before q-value filtering....

#df = df[df.fdr < fdr_treshold] # Filter Peptides with 

runs = df.Run.unique()

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


#eps = 1e-9
df_res = pd.concat(df_runs).reset_index()

df = df_res[df_res.fdr < fdr_treshold]

df_pivot = df.pivot(index="Stripped.Sequence", columns = "Run", values = "Precursor.Quantity")


