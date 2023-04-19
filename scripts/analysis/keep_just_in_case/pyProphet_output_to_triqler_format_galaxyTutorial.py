#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 00:26:04 2021

@author: ptruong
"""
import os 

import pandas as pd 
import numpy as np


# merge 
df = pd.DataFrame()
for i in os.listdir():
    if i[-4:] == ".tsv":
        if i[-7:] == "pqr.tsv":
            print("TEST")
        elif i[-7:] == "ins.tsv":
            print("TTT")
        else:
            df = df.append(pd.read_csv(i, sep = "\t"))

df.to_csv("merged_pyprophet.tsv", sep = "\t")
# map

def get_condition(x):
    condition_A = ["Sample"+str(i) for i in range(1,5)]
    condition_B = ["Sample"+str(i) for i in range(5,9)]
    
    cond = x.split("/")[1].split("_")[0]
    if cond in condition_A:
        return "A"
    elif cond in condition_B:
        return "B"
    else:
        return np.nan

condition_mapper = lambda x : get_condition(x)    

df_triq = df[["filename", "Charge", "Intensity", "m_score", "FullPeptideName", "ProteinName"]]

df_triq["condition"] = df.filename.map(condition_mapper)

df_triq = df_triq.rename(columns = {"filename":"run", "Charge":"charge", "Intensity":"intensity", 
                          "m_score":"searchScore", "FullPeptideName":"peptide", "ProteinName":"proteins"})
    
df_triq = df_triq[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]    
df_triq.searchScore = np.log(df_triq.searchScore)  # log-transform searchScores according to Matthiew  
df_triq.to_csv("triqler_format.csv", sep = "\t", index = False)
