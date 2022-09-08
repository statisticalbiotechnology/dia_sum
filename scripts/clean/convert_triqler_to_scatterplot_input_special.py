#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 20:05:17 2022

@author: ptruong
"""



import pandas as pd
import numpy as np
from parsers.parse_triqler import parse_triqler
import argparse

# params
iPipeline = "ID"
#iMethod = "msstats"
#iMethod = "msqrob2"
#iMethod = "top3"
#iMethod = "triqler"
triqler_results = f"results/{iPipeline}/triqler_results.csv"
method_results = f"results/{iPipeline}/{iMethod}_scatter_format_input.csv"


def convert_triqler_to_scatterplot_input_special(triqler_results, method_input):

    df_triqler_res = parse_triqler(triqler_results)
    df_triqler_res = (df_triqler_res[~((df_triqler_res.iloc[:,df_triqler_res.columns.str.contains("Pedro")] == 1).sum(axis=1) == 6)]) # Remove proteins that are not updates at any sample
    df_triqler_res.rename({"protein":"ProteinName"}, axis = 1, inplace = True)
    df_method_res = pd.read_csv(method_results, sep = "\t")
    df_method_res.set_index("ProteinName", inplace = True)
    # Find intersection of proteins
    
    intersecting_proteins = list(set(df_method_res.ProteinName) & set(df_triqler_res.ProteinName))
    
    df_triqler_res = df_triqler_res[df_triqler_res.ProteinName.isin(intersecting_proteins)]
    df_triqler_res["decoy"] = df_triqler_res.ProteinName.map(lambda x:x.split("_")[0])
    df_triqler_res = df_triqler_res[df_triqler_res["decoy"] != "DECOY"] #dilter away decoy
    df_triqler_res["specie"] = df_triqler_res.ProteinName.map(lambda x:x.split("_")[1])    
    protein_list = df_triqler_res.ProteinName.unique()
        
    df = pd.DataFrame(-df_triqler_res.set_index("ProteinName").log2_fold_change).rename({"log2_fold_change":"log2(A,B)"}, axis = 1)
    df["specie"] = df.index.map(lambda x:x.split("_")[1])
    df["q"] = df_triqler_res.set_index("ProteinName").q_value
    df["1"] = df_method_res["1"]
    df["2"] = df_method_res["2"]
    df.reset_index(inplace=True)    
    df = df[['specie', 'ProteinName', '1', '2', 'q', 'log2(A,B)']]

    return df

# ToDO
# Create a new snakemake rule for these specials
# Create argparse
# Use same plotting function and a special snakemake rule to create new plots


    
