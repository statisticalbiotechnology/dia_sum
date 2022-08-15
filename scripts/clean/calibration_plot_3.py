#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 15:08:28 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
from parsers.parse_triqler import  parse_triqler
import seaborn as sns
import matplotlib.pyplot as plt
import numpy.random as npr
import argparse
sns.set_context("talk")


  
triqler = parse_triqler(triqler_file)
triqler[~triqler.protein.str.contains("DECOY")]
triqler["specie"] = triqler.protein.map(lambda x:x.split("_")[1])
triqler["FDR"] = triqler["q_value"]
triqler.rename({"protein":"Protein"}, axis = 1, inplace = True)
top3 = pd.read_csv(top3_file, sep = "\t").rename({"q":"FDR", 'log2(A,B)':"log2FC"}, axis = 1)
top3.rename({"ProteinName":"Protein"}, axis = 1, inplace = True)
msstats = pd.read_csv(msstats_file, sep = ",").rename({"adj.pvalue":"FDR"}, axis = 1)
msstats["specie"] = msstats.Protein.map(lambda x:x.split("_")[1])
msqrob2.rename({"protein":"Protein"}, axis = 1, inplace = True)
msqrob2 = pd.read_csv(msqrob2_file, sep = ",").rename({"Unnamed: 0":"protein", "adjPval":"FDR", "logFC":"log2FC"}, axis = 1)
msqrob2["specie"] = msqrob2.protein.map(lambda x:x.split("_")[1])
msqrob2.rename({"protein":"Protein"}, axis = 1, inplace = True)

def get_fraction_hela(df_in, method):
    df=df_in.copy()
    df.sort_values(by = "FDR", inplace = True)
    df["count_HUMAN"] = df.Protein.str.contains("_HUMAN").astype(int)
    df["count_ECOLI"] = df.Protein.str.contains("_ECOLI").astype(int)
    df["count_YEAST"] = df.Protein.str.contains("_YEAST").astype(int)
    df["cumsum_HUMAN"] = df["count_HUMAN"].cumsum()
    df["cumsum_ECOLI"] = df["count_ECOLI"].cumsum()
    df["cumsum_YEAST"] = df["count_YEAST"].cumsum()
    df["Fraction_HeLa"] = df["cumsum_HUMAN"] / (df["cumsum_HUMAN"] + df["cumsum_ECOLI"] + df["cumsum_YEAST"])
    df["Method"] = method
    df.fillna(1, inplace = True) # assume NaN is maxed out FDR
    #df = df.reset_index().drop("index",axis = 1)
    return df[["FDR", "Fraction_HeLa", "Method"]]

def get_fraction_hela_df(triqler, top3, msstats, msqrob2):
    fdrs = []
    for df, method in zip([triqler, top3, msstats, msqrob2], ["Triqler", "Top3", "MSstats", "MSqRob2"]):
        fdrs.append(get_fraction_hela(df, method))
    df = pd.concat(fdrs)
    df = df.reset_index().drop("index", axis = 1).copy()
    return df








