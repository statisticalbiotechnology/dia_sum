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
import multiprocessing
from math import floor

os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811")

fdr_treshold = 0.01

df = pd.read_csv("report_recomputed_fdr.tsv", sep = "\t")
df.sort_values("CScore", ascending = False, inplace = True)
#df = df.groupby("Stripped.Sequence").head(1) # select top scoring PSM

# We need to recompute fdr before q-value filtering....

#df = df[df.fdr < fdr_treshold] # Filter Peptides with 

runs = df.Run.unique()

def compute_fdr(df_run, start_i, end_i, n_proc):
    df_decoy = df_run[df_run.decoy == True]
    df_target = df_run[df_run.decoy == False]
    cscores = []
    fdrs = []
    iteration = 0
    for i in df_run.index.unique()[start_i:end_i]:
        n_target = len(df_target[df_target.index > i])
        n_decoy = len(df_decoy[df_decoy.index > i])
        if (n_target + n_decoy) > 0:
            fdr_i = n_decoy/(n_target + n_decoy)
        fdrs.append(fdr_i)
        cscores.append(i)
        iteration += 1
        if iteration % 500 == 0:
            print(str(iteration) + " of " + str(len(df_run.index.unique())/n_proc))
    #fdr_map = dict(zip(cscores, fdrs))
    #df_run["fdr"] = df_run.index.map(fdr_map).fillna(0)
    return dict(zip(cscores, fdrs))

from time import time
start = time()
df_runs = []
for run in runs:
    df_run = df[df.Run == run].sort_values(by = "CScore")
    df_run.set_index("CScore", inplace = True)
    df_run = df_run.groupby("Stripped.Sequence").head(1) #select top scoring PSM
    #df_run = compute_fdr(df_run)
    df_runs.append(df_run)
    print(time()-start)
end = time()-start
print("Done " + str(end))

df_res = pd.concat(df_runs).reset_index()
df_res.set_index("CScore", inplace = True)

n_proc = 6
unique_index = len(df_res.index.unique())
step = floor(unique_index/n_proc)
residual = unique_index - n_proc*step
steps = [step*i for i in range(n_proc)]
steps.append(steps[-1]+step+residual)
compute_fdr_params = [(df_res, steps[i], steps[i+1], n_proc) for i in range(len(steps[:-1]))]

# This takes bout 20min for my diann dataset.
with multiprocessing.Pool() as pool:
    res = pool.starmap(compute_fdr, compute_fdr_params)


fdr_map = dict()
for i in range(len(res)):
    fdr_map.update(res[i])

df_run["fdr"] = df_run.index.map(fdr_map).fillna(0)






df = df_res[df_res.fdr < fdr_treshold]
df = df[df.decoy != True] # Filter away decoy as per msqrob2 manual
#df_pivot = df.pivot(index="Stripped.Sequence", columns = "Run", values = "Precursor.Quantity") #22892

df_pivot = df.pivot(index=["Stripped.Sequence","Protein.Ids"], columns = "Run", values = "Precursor.Quantity")
df_pivot.reset_index(inplace=True)
df_pivot["Proteins"] = df_pivot["Protein.Ids"]
df_pivot.drop("Protein.Ids", axis = 1)
df_pivot.to_csv("msqrob2_input_20220131.tsv", sep = "\t", index = False)


msq = pd.read_csv("msqrob2_input_20220131.tsv", sep = "\t")

# map sample

col = msq.columns[1:][1]

col_map = {}
for i in msq.columns:
    try:
        x = i.split("_")[8]
        print(x)
        if x == "1":
            new_col = "A_" + i
        if x == "2":
            new_col = "B_" + i
        col_map.update({i:new_col})
    except:
        col_map.update({i:i})

msq = msq.rename(columns=col_map)
msq.to_csv("msqrob2_input_20220131.tsv", sep = "\t", index = False)

p = pd.read_csv("msqrob2_input_20220131.tsv", sep = "\t")


# AFter running R code
test = pd.read_csv("msqrob2_results.tsv", sep = ",", index_col = 0)




