#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:23:55 2021

@author: ptruong
"""


import numpy as np
import pandas as pd

directory = "/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811"
filename = "report.tsv"
qvalue_treshold = 0.01
output_name = "triqler_input_diann_searchScore_Qvalue_treshold_0.01.csv"

def diann_to_triqler(filename, directory, output_name, qvalue_treshold = 1.00):
    df = pd.read_csv(directory + "/" + filename, sep = "\t")
    
    run_mapper = lambda x : x.split("_")[5]
    condition_mapper = lambda x : x.split("_")[8]
    
    df["run"] = df["Run"].map(run_mapper)
    df["condition"] = df["Run"].map(condition_mapper)
    df["charge"] = df["Precursor.Charge"]
    #df["searchScore"] = df["CScore"]
    df["searchScore"] = df["Q.Value"]
    df = df[df["searchScore"] < qvalue_treshold]
    df["searchScore"] = -np.log(df["searchScore"])
    df["intensity"] = df['Precursor.Quantity']
    df["peptide"] = df["Stripped.Sequence"]
    df["proteins"] = df["Protein.Ids"]
    
    df_triq = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    df_triq.to_csv(directory + "/" + output_name, sep = "\t", index = False)
    #return df_triq

if __name__ == "__main__":
    df_triq = diann_to_triqler(filename, directory, output_name, qvalue_treshold = qvalue_treshold) # fpr msqrobsum input
    
    # https://usermanual.wiki/Document/DIANN20GUI20manual.1528310561/view 
    
    
    