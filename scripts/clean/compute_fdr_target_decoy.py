#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:27:26 2022

@author: ptruong
"""

import pandas as pd 
import argparse
from math import floor
import multiprocessing

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

def target_decoy_fdr(input_file, output):
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
  df_res.to_csv(output, sep = "\t", index = False)
  
  
parser = argparse.ArgumentParser(
    description='Recomputes the FDR with target-decoy method.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input', type=str,
                    help='Input file to compute fdr.')

parser.add_argument('--output', type=str,
                    help='Output name.')

#parser.add_argument('--pcol', type=str,
#                    help='p-value or similar metric columns to be used to recompute the FDR.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input
output = args.output

if __name__ == "__main__":
    target_decoy_fdr(input_file, output)
    



  
