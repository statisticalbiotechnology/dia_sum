#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 09:20:02 2022

@author: ptruong
"""

import os
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

import pandas as pd
import numpy as np
import scipy.stats as stats
from q_value import qvalues
import multiprocessing
from math import floor

os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")


dfs = []
for file in os.listdir():
    if file[-10:] == "dscore.csv":
        dfs.append(pd.read_csv(file, sep = "\t"))
#        dfs.append(read_in_and_filter(file, m_score_treshold=0.01))
#        break
        #print(len(df_part))
        #df = pd.concat([df, df_part],axis = 1)        
df = pd.concat(dfs, axis = 0)
df.sort_values("m_score", ascending = True, inplace = True)

# For each run take top psm.
runs = df.filename.unique()


from time import time
start = time()
df_runs = []
for run in runs:
    df_run = df[df.filename == run].sort_values(by = "d_score")
    df_run.set_index("d_score", inplace = True)
    df_run = df_run.groupby("FullPeptideName").head(1) #select top scoring PSM
    #df_run = compute_fdr(df_run)
    df_runs.append(df_run)
    print(time()-start)
end = time()-start
print("Done " + str(end))

df_res = pd.concat(df_runs).reset_index()
df_res.set_index("d_score", inplace = True)
df_res.sort_values("d_score", ascending = False, inplace = True)
#df_peptide = df.groupby("FullPeptideName").head(1) #select top peptide as psm
#df_peptide = df_peptide.set_index("d_score") #d-score discriminate target from decoy
#df_peptide.sort_values("d_score", ascending = False, inplace = True) #reverse it again so that the logic of the code the same as the logic from the msqrob2-diann parser

def compute_fdr(df_run, start_i, end_i, n_proc):
    df_decoy = df_run[df_run.decoy == 1]
    df_target = df_run[df_run.decoy == 0]
    cscores = []
    fdrs = []
    iteration = 0
    for i in df_run.index.unique()[start_i:end_i]:
        n_target = len(df_target[df_target.index > i])
        n_decoy = len(df_decoy[df_decoy.index > i])
        if (n_target + n_decoy) > 0:
            fdr_i = n_decoy/(n_target + n_decoy)
        try:
            fdrs.append(fdr_i)
            cscores.append(i)
        except:
            pass
        iteration += 1
        if iteration % 500 == 0:
            print(str(iteration) + " of " + str(len(df_run.index.unique())/n_proc))
    #fdr_map = dict(zip(cscores, fdrs))
    #df_run["fdr"] = df_run.index.map(fdr_map).fillna(0)
    return dict(zip(cscores, fdrs))


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


df_res["fdr"] = df_res.index.map(fdr_map).fillna(0)

#df_peptide["fdr"] = df_peptide.index.map(fdr_map).fillna(0)
#df_res = df_peptide[df_peptide["fdr"] < 0.01]
#df_res = df_res[df_res.decoy == 0] #filter away decoy 

df_res = df_res[df_res["fdr"] < 0.01]
df_res = df_res[df_res.decoy == 0] #filter away decoy 

df_pivot = df_res.pivot(index=["FullPeptideName", "ProteinName"], columns = "filename", values = "Intensity")
df_pivot.reset_index(inplace=True)
df_pivot["Proteins"] = df_pivot["ProteinName"]
df_pivot.drop("ProteinName", axis = 1)
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



