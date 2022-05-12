#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:30:25 2022

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
        try:
            x = i.split("_")[8]
            #print(x)
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
    description='Converts diann to MsqRob2 input format. It also computes the FDR with target-decoy method.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name with recomputed q-values. Usually the report.tsv from DIANN')

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

