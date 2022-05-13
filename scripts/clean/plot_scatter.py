#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:22:40 2022

@author: ptruong
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from parsers.parse_triqler import parse_triqler



def convert_msqrob2_to_scatterplot(input_res, input_protein_quant):
    df = pd.read_csv(input_protein_quant, sep = ",").rename({'Unnamed: 0':"Proteins"}, axis = 1).set_index("Proteins")
    df_ = pd.read_csv(input_res, sep = ",").rename({'Unnamed: 0':"Proteins"}, axis = 1).set_index("Proteins")
    df = pd.concat([df ,df_], axis=1, join='inner')
    df["2"] = df.iloc[:,df.columns.str.contains("Sample_2")].mean(axis = 1)
    df["1"] = df.iloc[:,df.columns.str.contains("Sample_1")].mean(axis = 1)
    
    df = df[["1", "2", "pval", "adjPval", "logFC"]].reset_index()
    df["specie"] = df.Proteins.map(lambda x:x.split("_")[1])
    df = df[["specie", "Proteins", "1", "2", "pval", "adjPval", "logFC"]].rename({"Proteins": "ProteinName", 
                                                                                  "logFC":"log2(A,B)",
                                                                                  "pval":"p",
                                                                                  "adjPval":"q"}, axis = 1)
    return df

def convert_msstats_to_scatterplot(input_res, input_protein_quant):    
    df = pd.read_csv(input_protein_quant, sep = ",")
    df_ = pd.read_csv(input_res, sep = ",")
    df = df[["Protein", "LogIntensities", "GROUP"]]
    df = pd.concat([df_.set_index("Protein"), 
               pd.DataFrame(df[df.GROUP == 1].groupby("Protein").mean().LogIntensities).rename({"LogIntensities": "1"}, axis = 1),
               pd.DataFrame(df[df.GROUP == 2].groupby("Protein").mean().LogIntensities).rename({"LogIntensities": "2"}, axis = 1)],
              axis = 1)
    df.reset_index(inplace=True)
    df["specie"] = df.Protein.map(lambda x:x.split("_")[1])
    df = df[["specie", "Protein", "1", "2", "pvalue", "adj.pvalue", "log2FC"]].rename(
        {"Protein": "ProteinName", "log2FC":"log2(A,B)", "pvalue":"p", "adj.pvalue":"q"}, axis = 1)
    return df

def convert_triqler_to_scatterplot(input_file):
    df = parse_triqler(input_file)
    df["2"] = df.iloc[:,df.columns.str.contains("2:")].mean(axis = 1)
    df["1"] = df.iloc[:,df.columns.str.contains("1:")].mean(axis = 1)
    df["decoy"] = df.protein.map(lambda x:x.split("_")[0])
    df = df[df["decoy"] != "DECOY"] #dilter away decoy
    df["specie"] = df.protein.map(lambda x:x.split("_")[1])
    df = df.rename({"protein":"ProteinName", "q_value":"q", "log2_fold_change":"log2(A,B)"}, axis = 1)
    df = df[["specie", "ProteinName", "1", "2", "q", "log2(A,B)"]]
    df.specie.unique()
    return df



df_top3 = pd.read_csv("top3_results.csv", sep = "\t")


f, ax = plt.subplots(1, 1, figsize = (15,8))
#sns.scatterplot(ax = ax, data = df_final, x = "log2(B)", y = "log2(A,B)", alpha = 0.7, hue = "specie")
#df_final["specie"] = df_final.index.get_level_values("specie")
sns.regplot(data = df[df.specie == "ECOLI"], x = "1", y = "log2(A,B)", ax = ax, line_kws = {"ls":"--"}, label ="ECOLI")
sns.regplot(data = df[df.specie == "HUMAN"], x = "1", y = "log2(A,B)", ax = ax, line_kws = {"ls":"--"}, label = "HUMAN")
sns.regplot(data = df[df.specie == "YEAST"], x = "1", y = "log2(A,B)", ax = ax, line_kws = {"ls":"--"}, label = "YEAST")
ax.legend()
ax.axhline(1, linestyle = "--", color="tab:green", alpha = 0.5)
ax.axhline(0, linestyle = "--", color="tab:orange", alpha = 0.5)
ax.axhline(-2, linestyle = "--", color="tab:blue", alpha = 0.5)
#ax.grid()
#plt.legend(labels=['HUMAN (sample 1)', 'YEAST (sample 1)', "ECOLI (sample 1)", "HUMAN (sample 2)", "YEAST (sample 2)", "ECOLI (sample 2)"])
plt.title("std/mu ratio for log-transformed peptide values")


 
##### Converters



