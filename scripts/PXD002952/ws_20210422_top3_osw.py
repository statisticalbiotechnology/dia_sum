#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 21:39:55 2021

@author: ptruong
"""


import pandas as pd
import numpy as np

from scipy import stats

import os 

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/full_ts_v")

specie_mapper = lambda x : x.split("_")[-1]


def read_in_and_filter(filename, m_score_treshold = 0.01):  
    print(filename)
    df = pd.read_csv(filename, sep = "\t")
    df = df[df.decoy != 1]
    df = df[df.m_score < 0.01] # filter away crap, so all values should be good... we take average of top3 here
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    sample_id = df.sample_id[0]
    experiment_id = df.experiment_id[0]    
  
    def top3(df):
        df = (df.iloc[df['Intensity'].argsort()[::-1]]
                .groupby('ProteinName')['Intensity'].apply(lambda x: x.head(3).mean() if len(x.head(3)) >= 2 else np.nan)
                  .reset_index())
        #print(df.isna().sum())
        return df
    df_reduced = df[["ProteinName", "Intensity"]]
    df_protein = top3(df_reduced)
    df = df_protein
    df["specie"] = df.ProteinName.map(specie_mapper)
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df = df.set_index(["specie", "ProteinName"])
    df = pd.DataFrame(df.values, columns = midx, index = df.index)
    
    return df



# highest scoring vs top highest intenstiy peptide?
# We take top3 intensity because we filter away so much with m_score < 0.01, all values should be good.


df = pd.DataFrame()
for file in os.listdir():
    if file[-4:] == ".tsv":
        df_part = read_in_and_filter(file, m_score_treshold=0.01)
        print(len(df_part))
        df = pd.concat([df, df_part],axis = 1)        



df1 = df.iloc[:, df.columns.get_level_values("sample_id") == "1"]
df2 = df.iloc[:, df.columns.get_level_values("sample_id") == "2"]

(df1.isna().sum(axis = 1) < 3).sum() #approximately right amount of quantified peptied

(df2.isna().sum(axis = 1) < 3).sum() #approx right maount quantified peptides.


import time
import scipy.stats as stats


protein = df.index.get_level_values("ProteinName").unique()[0]

start = time.time()
protein_array = []
p_array  = []
proteins = df.index.get_level_values("ProteinName").unique()
for protein in proteins:
    df1_prots = df1.iloc[df1.index.get_level_values("ProteinName") == protein,:]
    df2_prots = df2.iloc[df2.index.get_level_values("ProteinName") == protein,:]
    p = stats.ttest_ind(df1_prots.T, df2_prots.T)[1]
    p_array.append(p[0])
    protein_array.append(protein)
end = time.time()

print(end-start)

# Significant dataframe.
df_p = pd.DataFrame(np.array([p_array, protein_array]).T, columns = ["p", "ProteinName"])
df_p["specie"] = df_p.ProteinName.map(specie_mapper)
df_p["p"] = df_p["p"].astype(float)
#df_p.index.name = "ProteinName"

# Dont need to join just find which proteins are significant and say index in the list on df

significant_proteins = df_p[df_p.p < 0.05].dropna()



#df.join(df_p, on = "ProteinName")
a = pd.DataFrame([1,2,3,4])
