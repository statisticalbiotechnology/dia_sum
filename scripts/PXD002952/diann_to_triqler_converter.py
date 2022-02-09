#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 23:08:38 2021

@author: ptruong
"""


import numpy as np
import pandas as pd


def diann_to_triqler(filename, qvalue_treshold = 1.00):
    df = pd.read_csv(filename, sep = "\t")
    
    run_mapper = lambda x : x.split("_")[5]
    condition_mapper = lambda x : x.split("_")[8]
    
    df["run"] = df["Run"].map(run_mapper)
    df["condition"] = df["Run"].map(condition_mapper)
    df["charge"] = df["Precursor.Charge"]
    #df["searchScore"] = df["CScore"]
    df["searchScore"] = df["Q.Value"]
    df["intensity"] = df['Precursor.Quantity']
    df["peptide"] = df["Stripped.Sequence"]
    df["proteins"] = df["Protein.Ids"]
    
    df_triq = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    return df_triq


if __name__ == "__main__":
    df_triq = diann_to_triqler("report.tsv")
    df_triq.intensity = np.log(df_triq.intensity) #if required to log intensity.
    df_triq.to_csv("triqler_input_diann.csv", sep = "\t", index = False)

