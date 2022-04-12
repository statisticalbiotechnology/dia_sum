#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 17:38:21 2022

@author: ptruong
"""

import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from q_value import qvalues
import multiprocessing
from math import floor
import argparse

def compute_fdr(df_run, start_i, end_i, n_proc):
    df_decoy = df_run[df_run.decoy == 1]
    df_target = df_run[df_run.decoy == 0]
    cscores = []
    fdrs = []
    iteration = 0
    for i in df_run.index.unique()[start_i:end_i]:
        n_target = len(df_target[df_target.index < i])
        n_decoy = len(df_decoy[df_decoy.index < i])
        if (n_target + n_decoy) > 0:
            fdr_i = n_decoy/(n_target)
        try:
            fdrs.append(fdr_i)
            cscores.append(i)
        except:
            pass
        iteration += 1
        if iteration % 500 == 0:
            print(str(iteration) + " of " + str(len(df_run.index.unique())/n_proc))
    return dict(zip(cscores, fdrs))

def convert_osw_to_msqrob2(input_file, output, fdr_threshold = 0.01):
    df = pd.read_csv(input_file, sep = "\t")
    
    # For each run take top psm.
    runs = df.filename.unique()
    
    from time import time
    start = time()
    df_runs = []
    for run in runs:
        df_run = df[df.filename == run].sort_values(by = "m_score")
        df_run.set_index("m_score", inplace = True)
        df_run = df_run.groupby("FullPeptideName").head(1) #select top scoring PSM
        df_runs.append(df_run)
        print(time()-start)
    end = time()-start
    print("Done " + str(end))
    
    df_res = pd.concat(df_runs).reset_index()
    df_res.set_index("m_score", inplace = True)
    df_res.sort_values("m_score", ascending = True, inplace = True)
        
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
    
    df_res["fdr"] = df_res.index.map(fdr_map).fillna(0)
    df_res = df_res[df_res["fdr"] < fdr_threshold]
    df_res = df_res[df_res.decoy == 0] #filter away decoy 
    
    df_pivot = df_res.pivot(index=["FullPeptideName", "ProteinName"], columns = "filename", values = "Intensity")
    df_pivot.reset_index(inplace=True)
    df_pivot["Proteins"] = df_pivot["ProteinName"]
    df_pivot.drop("ProteinName", axis = 1)
    msq = df_pivot
    
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
    msq.to_csv(output, sep = "\t", index = False)
        
 
parser = argparse.ArgumentParser(
    description='Convert OSW to MsqRob2 input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--fdr_threshold', type=float, default = 0.01,
                    help='fdr threshold to apply on qCol.')

parser.add_argument('--output', type=str, default = "osw_msqrob2_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
fdr_threshold = args.fdr_threshold
output = args.output
      
if __name__ == "__main__":
    convert_osw_to_msqrob2(input_file, output, fdr_threshold)
    
