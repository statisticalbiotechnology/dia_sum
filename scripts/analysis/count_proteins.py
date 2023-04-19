#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 13:28:18 2022

@author: ptruong
"""

import os
import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/dia_sum/scripts/clean/data/ID/osw_results")

for i in os.listdir():
    print(i)
    df = pd.read_csv(i, sep = "\t", usecols = ["decoy", "ProteinName", "FullPeptideName", "filename", "m_score"])
    df = df[df.decoy != 1]
    #df = df[df.m_score < 0.000794328234724281] # m_score cutoff with target FDR 0.01
    n_proteins = len(df.ProteinName.unique())
    n_peptides = len(df.FullPeptideName.unique())
    filename = df.filename.unique()
    print(filename)
    print("proteins:" + str(n_proteins))
    print("peptides:" + str(n_peptides))
    print("----")

os.chdir("/home/ptruong/git/dia_sum/scripts/clean/data/PS")

df = pd.read_csv("report.tsv", sep = "\t")
df.columns
df["Protein.Ids"] = df["Protein.Ids"].astype(str)
df = df[df["Protein.Ids"].map(lambda x:x.split("_")[0]) != "DECOY"]

for i in df["File.Name"].unique():
    df_i = df[df["File.Name"] == i]
    print(i)
    print("proteins: " + str(len(df_i["Protein.Ids"].unique())))
    print("peptides: " +str(len(df_i["Stripped.Sequence"].unique())))
    print("----")

df["PG.Q.Value"].max()













