#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:26:00 2022

@author: ptruong
"""

import os
import pandas as pd
from parsers.parse_triqler import  parse_triqler


file_dir = "/hdd_14T/data/PXD002952/20210614_dataset/result_files_20220214/PS/"

negCol = "Differential HeLa"
posCol = "Differential non-HeLa"

def remove_decoy(df, protein_column):   
    res = df[~df[protein_column].str.contains("DECOY_")]
    return res

def countProteins(df, method):
    negCol = "Differential HeLa"
    posCol = "Differential non-HeLa"
    df = remove_decoy(df, protein_column = "Protein")
    df.dropna(subset=["FDR"], inplace = True)
    df.sort_values(by = "FDR", inplace = True)
    negative = df["Protein"].str.contains("_HUMAN").astype(int)
    df[negCol] = negative.cumsum()
    df[posCol] = (1-negative).cumsum()
    df["method"] = method
    return df


triqler_results = parse_triqler("triqler_results/fc_0.92")
top3_results = pd.read_csv("top3_output_diann.csv", sep = "\t")
msstats_results = pd.read_csv("msstat_output.csv", sep = ",")
msqrob2_results = pd.read_csv("msqrob2_results.tsv", sep = ",").rename({"Unnamed: 0":"Protein"},axis=1)

#Rename protein, fdr to same name
triqler_results = triqler_results.rename({"q_value":"FDR", "protein":"Protein"}, axis = 1)
top3_results = top3_results.rename({"q":"FDR", "ProteinName":"Protein"}, axis = 1)
msstats_results = msstats_results.rename({"adj.pvalue":"FDR"}, axis = 1)
msqrob2_results = msqrob2_results.rename({"adjPval":"FDR"}, axis = 1)



methods = ["Triqler", "Top3", "MsStats", "MsqRob2"]
data = [triqler_results, top3_results, msstats_results, msqrob2_results]
zipped = zip(methods, data)


#map different protein names

dfs = []
for method, df in zipped:
    df_count = countProteins(df, method)
    dfs.append(df_count.loc[:,["Protein", posCol, negCol, "method"]])



dfs[0].columns
pd.concat(dfs)
    
df

plt.plot(df[negCol], df[posCol])

triqler_results

def countProteins(df, pCol, methodName):
    df.drop(df[df[pCol].str.contains("DECOY_")].index, inplace=True)
    negative = df[pCol].str.contains("_HUMAN").astype(int)
    df[negCol] = negative.cumsum()
    df[posCol] = (1-negative).cumsum()
    df[eCol] = df[negCol]/(df[posCol]+df[negCol])
    df["method"] = methodName
    return df








