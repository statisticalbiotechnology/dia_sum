#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:09:47 2021

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

    #midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df["specie"] = df["ProteinName"].map(specie_mapper)
    df = df.set_index(["specie", "ProteinName", "FullPeptideName", "sample_id", "experiment_id"])
    df = df[["Intensity"]]
    print(sample_id)
    print(experiment_id)
    df = pd.DataFrame(df.values, index = df.index)
    
    return df


dfs = []
for file in os.listdir():
    if file[-4:] == ".tsv":
        dfs.append(read_in_and_filter(file, m_score_treshold=0.01))
        #print(len(df_part))
        #df = pd.concat([df, df_part],axis = 1)        
df = pd.concat(dfs, axis = 0)
df = np.log2(df)

# WE CAN SPLIT FOR SAMPLE on df
# We CAN SPLIT FOR SPECIE

import matplotlib.pyplot as plt
import seaborn as sns
len(list(df.index.get_level_values("FullPeptideName").unique()))

def get_peptide_mu_sigma(df):
    df_means = df.groupby(df.index.get_level_values("FullPeptideName")).mean()
    
    df_stat = pd.DataFrame(df_means.values, index = df_means.index, columns = ["mu"])
    df_stat["std"] = df.groupby(df.index.get_level_values("FullPeptideName")).std()
    df_stat["std/mu-ratio"] = df_stat["std"] / df_stat["mu"]
    return df_stat




df_stat = get_peptide_mu_sigma(df)
f, ax = plt.subplots(1, 1, figsize = (15,8))
sns.scatterplot(ax = ax, data = df_stat, x = "mu", y = "std/mu-ratio", alpha = 0.5)
plt.title("std/mu ratio for log-transformed peptide values")

def select_specie_and_sample(df, specie = "HUMAN", sample = "1"):
    return df.iloc[(df.index.get_level_values("specie") == specie) & (df.index.get_level_values("sample_id") == sample), :]


get_peptide_mu_sigma(select_specie_and_sample(df, specie = "HUMAN", sample = "1"))
get_peptide_mu_sigma(select_specie_and_sample(df, specie = "YEAS8", sample = "1"))
get_peptide_mu_sigma(select_specie_and_sample(df, specie = "ECOLI", sample = "1"))
get_peptide_mu_sigma(select_specie_and_sample(df, specie = "HUMAN", sample = "2"))
get_peptide_mu_sigma(select_specie_and_sample(df, specie = "YEAS8", sample = "2"))
get_peptide_mu_sigma(select_specie_and_sample(df, specie = "ECOLI", sample = "2"))

f, ax = plt.subplots(1, 1, figsize = (15,8))
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "HUMAN", sample = "1")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "YEAS8", sample = "1")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "ECOLI", sample = "1")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "HUMAN", sample = "2")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "YEAS8", sample = "2")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
sns.scatterplot(ax = ax, data = get_peptide_mu_sigma(select_specie_and_sample(df, specie = "ECOLI", sample = "2")), x = "mu", y = "std/mu-ratio", alpha = 0.3)
plt.legend(labels=['HUMAN_1', 'YEAST_1', "ECOLI_1", "HUMAN_2", "YEAST_2", "ECOLI_2"])
plt.title("std/mu ratio for log-transformed peptide values")


# pandas count
sns.histplot(df.groupby(df.index.get_level_values("FullPeptideName")).count())

df_counts = []
for specie in ["HUMAN", "YEAS8", "ECOLI"]:
    for sample in ["1", "2"]:
        df_tmp = select_specie_and_sample(df, specie = specie, sample = sample)
        df_tmp = df_tmp.groupby(df_tmp.index.get_level_values("FullPeptideName")).count() 
        df_tmp = pd.DataFrame(df_tmp.values, index = df_tmp.index, columns = ["count"])
        df_tmp["specie"] = specie
        df_tmp["sample"] = sample
        df_tmp["specie-sample"] = specie + "_" + sample
        df_counts.append(df_tmp)
df_count = pd.concat(df_counts, axis = 0).reset_index()

ax = sns.countplot(x="count", hue="specie-sample", data=df_count)




