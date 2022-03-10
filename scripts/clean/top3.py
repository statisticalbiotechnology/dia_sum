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

#from parsers.parse_triqler import parse_triqler
from q_value import qvalues
#from mappers import experiment_id_mapper, sample_id_mapper, specie_mapper, decoy_mapper
#import top3_formatter.diann
#import top3_formatter.osw

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





# DIANN
df = pd.read_csv("report_for_top3.tsv", sep = "\t")
#triq = pd.read_csv("triqler_input_diann_searchScore_Qvalue.csv", sep = "\t")                      
fdr_treshold = 0.01
df = df[df["Q.Value"] < fdr_treshold]
df = top3_formatter.diann.avg_3_largest_precursor(df)
df = top3(df)


# OSW
df = top3_formatter.osw.avg_3_largest_precursor("osw_results", m_score_treshold = 0.01)
df = top3(df)




















