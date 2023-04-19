#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 19:29:52 2021

@author: ptruong
"""

import pandas as pd
import os 
os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")

dfs = []
for file in os.listdir():
    if file[-10:] == "dscore.csv":
        df = pd.read_csv(file, sep = "\t")
        df = df[df["m_score"] < 0.01]
        dfs.append(df)

df_concat = pd.concat(dfs)
df_concat.to_csv("concatenated_osw_output_m_score_0.01.csv", sep = "\t")








