#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:23:43 2021

@author: ptruong
"""


import os 
import pandas as pd

def generate_nfcore_input(filepath):

    files = []
    for file in os.listdir(filepath):
        if file[-5:] == ".mzML":
            files.append(file)
            
    df = pd.DataFrame(files, columns = ["Spectra_Filepath"])  
    
    x = df.iloc[0][0]
    
    sample_mapper = lambda x: x.split("_")[5]
    condition_mapper = lambda x: int(x.split("_")[8])
    bioReplicate_mapper = lambda x: int(x.split("_")[-1].split(".")[0][-1])
    
    df = df.reset_index().rename({"index": "Sample"}, axis = "columns")
    df["Sample"] += 1
    df["BatchID"] = "LFQ"
    df["MSstats_Condition"] = df.Spectra_Filepath.map(condition_mapper)
    df["MSstats_BioReplicate"] = df.Spectra_Filepath.map(bioReplicate_mapper)
    df[["Sample", "BatchID", "MSstats_Condition", "MSstats_BioReplicate", "Spectra_Filepath"]]
    return df

def generate_nfcore_spectral_library_input():
    headers = ["Sample" , "BatchID", "Library_Filepath"]
    data = ["1", "LFQ", "library.tsv"]
    df = pd.DataFrame(data, index = headers).T
    return df

def generate_nfcore_irts_input():
    headers = ["Sample" , "BatchID", "Library_Filepath"]
    data = ["1", "LFQ", "hroest_DIA_iRT.TraML"]
    df = pd.DataFrame(data, index = headers).T
    return df

df = generate_nfcore_input(os.getcwd())
df.to_csv("sample_input.csv", sep = "\t", index = False)

lib_df = generate_nfcore_spectral_library_input()
lib_df.to_csv("input_spectral_library.csv", sep = "\t", index = False)

generate_nfcore_irts_input()


