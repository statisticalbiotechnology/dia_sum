#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:19:23 2022

@author: ptruong
"""


import os 
import pandas as pd
import numpy as np
import scipy.stats as stats

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from parsers.parse_triqler import parse_triqler
from q_value import qvalues
from mappers import experiment_id_mapper, sample_id_mapper, specie_mapper, decoy_mapper


# this is just wierd?
def compute_triqler_top3_submodule(run, fdr_treshold = 0.01):
    triq_run = triq[triq.run == run]   
    triq_run.searchScore = np.e**(-1 * triq_run.searchScore) #re-convert searchScore, because our triqler_input has -np.log(-Q.Value) from diann result.tsv as searchScore
    triq_run = triq_run[triq_run["searchScore"] < fdr_treshold]
    
    def triqler_top3(triq_run):
        res = triq_run.groupby("proteins")["intensity"].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan).reset_index()
        return res
    
    def triqler_printout_unique_peptides_proteins(run):
        condition = triq_run.condition.unique()
        print(f"run : {run} - condition : {condition}")
        print(f"Unique peptides detected: {len(triq_run.peptide.unique())}")
        print(f"Unique proteins detected: {len(triq_run.proteins.unique())}")
        print()
    
    triqler_printout_unique_peptides_proteins((run))
    res = triqler_top3(triq_run)
    experiment_id = triq_run.run.unique()[0]
    sample_id = triq_run.condition.unique()[0]
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df = pd.DataFrame(res.values)  
    specie_map = lambda x: x.split("_")[-1]
    protein_map = lambda x: x.split("_")[-2]
    res["specie"] = res.proteins.map(specie_map)
    res["ProteinName"] = res.proteins
    res = res.drop(["proteins"], axis = 1)
    res = res.set_index(["specie", "ProteinName"])
    df = pd.DataFrame(res.values, columns = midx, index = res.index)
    return df

def avg_3_largest_precursor_on_run_level(df_run):
    def avg_of_3_largest_precursor(df_run):
        return df_run.groupby("Protein.Ids")["Precursor.Quantity"].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan).reset_index()
    
    res = avg_of_3_largest_precursor(df_run)
    experiment_id = df_run.Run.unique()[0].split("_")[5]
    sample_id = df_run.Run.unique()[0].split("_")[8]
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df = pd.DataFrame(res.values)
    specie_map = lambda x: x.split("_")[-1]
    protein_map = lambda x: x.split("_")[-2]
    res["specie"] = res["Protein.Ids"].map(specie_map)
    res["ProteinName"] = res["Protein.Ids"]
    res = res.drop(["Protein.Ids"], axis = 1)
    res = res.set_index(["specie", "ProteinName"])
    df = pd.DataFrame(res.values, columns = midx, index = res.index)
    return df

def avg_3_largest_precursor(df):
    dfs = []
    for run in df.Run.unique(): # for every tun otherwise we mix together proteins from multiple runs when doing top3.
        df_run = df[df.Run == run]
        df_top3_run_level = avg_3_largest_precursor_on_run_level(df_run)
        dfs.append(df_top3_run_level)
    df = pd.concat(dfs, axis = 1)
    return df

def top3(df):
    A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "1"].isna().sum(axis=1)<2]
    A = A.iloc[:,A.columns.get_level_values("sample_id") == "1"]
    B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "2"].isna().sum(axis=1)<2]
    B = B.iloc[:,B.columns.get_level_values("sample_id") == "2"]
    
    # Find overlapping proteins
    overlapping_proteins = list(set(A.index) & set(B.index))
    A = A[A.index.isin(overlapping_proteins)]
    B = B[B.index.isin(overlapping_proteins)]
    
    p_vals = stats.ttest_ind(A, B, axis = 1)[1]
    p_vals = pd.DataFrame(p_vals, columns = ["p"])
    p_vals["q"] = qvalues(p_vals)
    p_vals = pd.DataFrame(p_vals.values, index = A.index, columns = ["p", "q"])
    p_vals.sort_values("q",inplace =True)
    p_vals = p_vals.astype(float)
    
    A = np.log2(A.sum(axis=1))
    B = np.log2(B.sum(axis=1))
    
    A.name = "1"
    B.name = "2"
    
    df_final = pd.concat([A, B, p_vals], axis = 1)
    df_final["log2(A,B)"] = df_final["1"] - df_final["2"]
    return df_final

df = pd.read_csv("report.tsv", sep = "\t")
#triq = pd.read_csv("triqler_input_diann_searchScore_Qvalue.csv", sep = "\t")                      
fdr_treshold = 0.01
df = df[df["Q.Value"] < fdr_treshold]
df = avg_3_largest_precursor(df)
df = top3(df)

