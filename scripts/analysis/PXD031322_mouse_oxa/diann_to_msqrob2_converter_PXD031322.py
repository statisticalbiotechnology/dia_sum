#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 12:04:12 2022

@author: ptruong
"""

import os
import pandas as pd
import numpy as np
import multiprocessing
from math import floor
import argparse

input_file = "report.tsv"
output = "msqrob2_input_20220409.tsv"
fdr_treshold = 0.01


def compute_fdr(df_run, start_i, end_i, n_proc):
    df_decoy = df_run[df_run.decoy == True]
    df_target = df_run[df_run.decoy == False]
    cscores = []
    fdrs = []
    iteration = 0
    for i in df_run.index.unique()[start_i:end_i]:
        n_target = len(df_target[df_target.index > i])
        n_decoy = len(df_decoy[df_decoy.index > i])
        if (n_target ) > 0:
            fdr_i = n_decoy/(n_target)
            fdrs.append(fdr_i)
            cscores.append(i)
        else:
            fdrs.append(1)
            cscores.append(i)
        iteration += 1
        if iteration % 500 == 0:
            print(str(iteration) + " of " + str(len(df_run.index.unique())/n_proc))
    #fdr_map = dict(zip(cscores, fdrs))
    #df_run["fdr"] = df_run.index.map(fdr_map).fillna(0)
    return dict(zip(cscores, fdrs))

def convert_diann_to_msqrob2(input_file, output, fdr_threshold):
    """
    Converts and computes the FDR using target-decoy method.
    """
    df = pd.read_csv(input_file, sep = "\t")
    
    # map float in Protein.Ids to check what is the problem and remove these
    float_mapper = lambda x : isinstance(x, float)
    df["isfloat"] = df["Protein.Ids"].map(float_mapper)
    df[df["isfloat"] == True]["Protein.Ids"]
    df = df[df["isfloat"] != True]

    # get decoy column
    decoy_mapper = lambda x:x.split("_")[0] == "DECOY"
    df["decoy"] = df["Protein.Ids"].map(decoy_mapper)
    
    df.sort_values("CScore", ascending = False, inplace = True)
    runs = df.Run.unique()
    
    
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
    df_res.sort_values("CScore", ascending = False, inplace = True)
    
    n_proc = 7
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
    df_res["q"] = df_res.index.map(fdr_map).fillna(0)
    df = df_res[df_res["q"] < fdr_treshold]
    df = df[df.decoy != True] # Filter away decoy as per msqrob2 manual
    #df_pivot = df.pivot(index="Stripped.Sequence", columns = "Run", values = "Precursor.Quantity") #22892
    df_pivot = df.pivot(index=["Stripped.Sequence","Protein.Ids"], columns = "Run", values = "Precursor.Quantity")
    df_pivot.reset_index(inplace=True)
    df_pivot["Proteins"] = df_pivot["Protein.Ids"]
    df_pivot.drop("Protein.Ids", axis = 1)
    
    msq = df_pivot
    
    # map sample    
    col_map = {}
    for i in msq.columns:
        #print(i)
        try:
            x = i.split("-")[-1][:-1]
            #print(x + " ---- " + i)
            if x == "ST":
                new_col = "ST_" + i
                col_map.update({i:new_col})
            elif x == "LT":
                new_col = "LT_" + i
                col_map.update({i:new_col})
            elif x == "Ctrl":
                new_col = "Ctrl_" + i
                col_map.update({i:new_col})
            else:
                col_map.update({i:i})
        except:
            col_map.update({i:i})
    #col_map["proteins"] = "protein"
    msq = msq.rename(columns=col_map)
    msq = msq.T[~msq.columns.str.contains("GPF")].T # drop GPF columns
    
    # All msqrob2 examples in vignette är case control examples so i split up the dataset.
    msq.T[~msq.columns.str.contains("ST")].T.to_csv("Ctrl_LT_"+output, sep = "\t", index = False) #Ctrl-LT
    msq.T[~msq.columns.str.contains("LT")].T.to_csv("Ctrl_ST_"+output, sep = "\t", index = False) #Ctrl-ST
    msq.T[~msq.columns.str.contains("Ctrl")].T.to_csv("LT_ST_"+output, sep = "\t", index = False) #LT-ST
    msq.to_csv("LT_ST_Ctrl_" + output, sep = "\t", index = False)


parser = argparse.ArgumentParser(
    description='Converts diann to MsqRob2 input format. It also computes the FDR with target-decoy method.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name. Usually the report.tsv from DIANN')

parser.add_argument('--fdr_threshold', type=float, default = 0.01,
                    help='fdr threshold to apply on q-value.')

parser.add_argument('--output', type=str, default = "msqrob2_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
fdr_threshold = args.fdr_threshold
output = args.output

if __name__ == "__main__":
    convert_diann_to_msqrob2(input_file, output, fdr_threshold)

