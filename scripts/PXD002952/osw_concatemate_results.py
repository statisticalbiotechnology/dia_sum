#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 06:51:20 2022

@author: ptruong
"""

import os 
os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")

import pandas as pd
from time import time


start = time()
dfs = []
for i in os.listdir():
    if i[-10:] == "dscore.csv":
        df = pd.read_csv(i, sep = "\t")
        dfs.append(df)
        print(time()-start)
end = time()
print(end-start)

df = pd.concat(dfs)
df.to_csv("concatenated_osw_results.csv", sep = "\t")

# Compute m_score cut-off with Swath2Stats
# m_score cutoff

m_score_cutoff = 0.00079433
df_tresh = df[df.m_score < m_score_cutoff]

df_example = pd.read_csv("example_disaggregate_output_format.csv", sep = ",")




# count transitions

df.reset_index(inplace=True)
x = df["aggr_Fragment_Annotation"][0]

count_transitions = lambda x:len(x.split(";"))


df["transition_count"] = df["aggr_Fragment_Annotation"].map(count_transitions)

df["transition_count"].max() #92
df["transition_count"].min() #4


df.transition_count
df.columns

df["masserror_ppm"]


transition_group_id_count = df.groupby("transition_group_id").count()

transition_group_id_count.transition_count.min()







