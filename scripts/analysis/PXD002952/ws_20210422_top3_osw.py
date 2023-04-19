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

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

from q_value import qvalues
from triqler_output_to_df import  parse_triqler
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/full_ts_v")


# filename has different formatting, we need to change number or implement regex.
experiment_id_mapper = lambda x: x.split("_")[5]
sample_id_mapper = lambda x: x.split("_")[8] #hye124 
specie_mapper = lambda x: x.split("_")[-1]

def read_in_and_filter(filename, m_score_treshold = 0.01):  
    print(filename)
    df = pd.read_csv(filename, sep = "\t")
    df = df[df.decoy != 1]
    df = df[df.m_score < m_score_treshold] # filter away crap, so all values should be good... we take average of top3 here
    print(str(len(df)) + " significantly identified peptides at " + str(m_score_treshold) + " FDR-treshold.")
    print("")
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    sample_id = df.sample_id[0]
    experiment_id = df.experiment_id[0]     
    def top3(df):
        df = (df.groupby('ProteinName')['Intensity'].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan)
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

dfs = []
for file in os.listdir():
    if file[-4:] == ".tsv":
        dfs.append(read_in_and_filter(file, m_score_treshold=0.01))
        #print(len(df_part))
        #df = pd.concat([df, df_part],axis = 1)        
df = pd.concat(dfs, axis = 1)


A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "1"].isna().sum(axis=1)<2]
A = A.iloc[:,A.columns.get_level_values("sample_id") == "1"]
#A = np.log2(A)
B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "2"].isna().sum(axis=1)<2]
B = B.iloc[:,B.columns.get_level_values("sample_id") == "2"]
#B = np.log2(B)

# Check the histograms
A.stack().hist(bins = 100)
B.stack().hist(bins = 100)

# Find overlapping proteins
overlapping_proteins = list(set(A.index) & set(B.index))
A = A[A.index.isin(overlapping_proteins)]
B = B[B.index.isin(overlapping_proteins)]

import scipy.stats as stats

p_vals = stats.ttest_ind(A, B, axis = 1)[1]
p_vals = pd.DataFrame(p_vals, columns = ["p"])
p_vals["q"] = qvalues(p_vals)
p_vals = pd.DataFrame(p_vals.values, index = A.index, columns = ["p", "q"])
p_vals.sort_values("q",inplace =True)
p_vals = p_vals.astype(float)

sns.lineplot(x = "p", y = "q", data=p_vals)



def plot_n_DE(p_vals):
    n_array = []
    for q in np.arange(0,1.01,0.01):
        n = (p_vals.q < q).sum()
        n_array.append(n)
    
    n_de = pd.DataFrame([n_array, np.arange(0,1.01,0.01)], index = ["n", "q_treshold"])
    plt.plot(np.arange(0,1.01,0.01), n_array)

#plot_n_DE(p_vals[p_vals.index.get_level_values("specie")=="HUMAN"])
#A = A.sum(axis=1)
#B = B.sum(axis=1)


A = np.log2(A.sum(axis=1))
B = np.log2(B.sum(axis=1))

A.name = "1"
B.name = "2"

df_final = pd.concat([A, B, p_vals], axis = 1)



df_final["log2(A,B)"] = df_final["1"] - df_final["2"]
df_h = df_final.iloc[df_final.index.get_level_values("specie") == "HUMAN", :]
df_y = df_final.iloc[df_final.index.get_level_values("specie") == "YEAS8", :]
df_e = df_final.iloc[df_final.index.get_level_values("specie") == "ECOLI", :]


import matplotlib.pyplot as plt
import seaborn as sns
f, ax = plt.subplots(1, 1)
sns.scatterplot(ax = ax, data = df_h, x = "2", y = "log2(A,B)", alpha = 0.5)
sns.scatterplot(ax = ax, data = df_y, x = "2", y = "log2(A,B)", alpha = 0.5)
sns.scatterplot(ax = ax, data = df_e, x = "2", y = "log2(A,B)", alpha = 0.5)
plt.legend(labels=['HUMAN', 'YEAS8', 'ECOLI'])








def plot_n_DE(df):
    n_array = []
    for fc in np.arange(0,2.05,0.05):
        n_gt = (df["log2(A,B)"] > fc).sum()
        n_lt = (df["log2(A,B)"] < -fc).sum()
        n = n_gt + n_lt
        df["log2(A,B)"].between(-fc, fc, inclusive=True)
        n_array.append(n)
    
    n_de = pd.DataFrame([n_array, np.arange(0,2.05,0.05)], index = ["n", "fc_treshold"])
    plt.plot(np.arange(0,2.05,0.05), n_array)

df_final_q = df_final[df_final["q"]<0.01]

plot_n_DE(df_final_q.iloc[df_final_q.index.get_level_values("specie") == "HUMAN",:])
plot_n_DE(df_final_q.iloc[df_final_q.index.get_level_values("specie") == "YEAS8",:])
plot_n_DE(df_final_q.iloc[df_final_q.index.get_level_values("specie") == "ECOLI",:])
plt.legend(labels=['HUMAN', 'YEAS8', 'ECOLI'])



os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_results")

fc_tresh = []
n_hs = []
n_ye = []
n_ec = []

for file in sorted(os.listdir()):       
    fc = float(file.split("_")[1])
    df_triq = parse_triqler(file)
    df_triq = df_triq[df_triq.q_value < 0.01]
    df_triq["specie"] = df_triq.protein.map(specie_mapper)
    n_hs.append((df_triq["specie"] == "HUMAN").sum())
    n_ye.append((df_triq["specie"] == "YEAS8").sum())
    n_ec.append((df_triq["specie"] == "ECOLI").sum())
    fc_tresh.append(fc)


plt.plot(fc_tresh, n_hs)
plt.plot(fc_tresh, n_ye)
plt.plot(fc_tresh, n_ec)
plt.legend(labels=['HUMAN', 'YEAS8', 'ECOLI'])

