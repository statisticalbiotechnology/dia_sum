#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 15:57:52 2021

@author: ptruong
"""

import os 
import pandas as pd


df_merged = pd.DataFrame()
for i in os.listdir():
    if i[-4:] == ".tsv":
        parsed = i.split("_")
        experiment_id = parsed[5]
        sample_id = parsed[7]
        df = pd.read_csv(i, sep = "\t")
        df["experiment_id"] = experiment_id
        df["sample_id"] = sample_id
        df_merged = df_merged.append(df)


df_merged.columns
df_merged
df_t = df_merged[df_merged.m_score < 0.01] #FDR-treshold

df = df_t[df_t.decoy == 0]

df.columns
df.experiment_id.unique()


#top 3 peptides to select protein
df_sample = df[df.experiment_id == "001-Pedro"]

protein = df_sample.ProteinName.unique()[1]

df_prot = df_sample[df_sample.ProteinName == protein]

if len(df_prot) > 3:
    # Take the sum of this as protein quantity
else:
    # if less than 3 peptides for protein... discard it.








































