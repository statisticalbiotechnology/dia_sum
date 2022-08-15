#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 18:49:20 2022

@author: ptruong
"""

# Script for computing top3, p-value and q-value for COVID-19 dataset



import pandas as pd
import numpy as np
import os 
import scipy.stats as stats

os.chdir("/home/ptruong/git/dia_sum/scripts/clean")
from q_value import qvalues

os.chdir("/hdd_14T/data/pxd021197/ftp.ebi.ac.uk/pride-archive/2021/07/PXD021197/DIA")

dfs = []
for i in os.listdir():
    if i[-4:] == ".tsv":
        if i[:2] == "CK":
            df = pd.read_csv(i, sep = "\t")
            dfs.append(df)

df = pd.concat(dfs)
df.to_csv("result.tsv", sep = "\t", index = False)


# Add sample_id and experiment_id
os.chdir("/hdd_14T/data/pxd021197")
mapper_file = "PXD021197_DIA_DDA.xlsx"
mapper = pd.read_excel(mapper_file)

os.chdir("/hdd_14T/data/pxd021197/ftp.ebi.ac.uk/pride-archive/2021/07/PXD021197/DIA")

filename = "result.tsv"
df = pd.read_csv(filename, sep = "\t", usecols = ["Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName", "filename"])
df["mapper"] = df.filename.map(lambda x:x.split(".")[0])
mapper.SampleName
len(mapper.Patient.unique())
mapper.columns
patient_map=mapper[["FileName", "Patient"]].set_index("FileName").to_dict()["Patient"]
run_map=mapper[["FileName", "SampleName"]].set_index("FileName").to_dict()["SampleName"]
df["patient"] = df.mapper.map(patient_map)
#df.condition = df.condition.map(control_condition_rename)
# case control rename
#    df.condition = df.condition.map(control_case_rename)
df["run"] = df.mapper.map(run_map)
#df["condition"] = df.run.map(lambda x:x.split("_")[2].split("+")[0].split("-")[0])
df["condition"] = df.run.map(lambda x:x.split("_")[2]) # every day split
df["condition"] = df.condition.map(lambda x:x.split("+")[0])
df["condition"] = df.condition.map(lambda x:x.split("-")[0])
df["experiment_id"] = df["run"]
df["sample_id"] = df["condition"]
df = df[df["sample_id"] != "NA"]

df["decoy"] = df.ProteinName.str.contains("DECOY")
df = df[df.experiment_id.str.contains("Plasma")]
# Top3 formatter
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/top3_formatter")



def avg_3_largest_precursor_on_run_level(df):  
    df = df[df.decoy != 1]
    #df = df[df.m_score < m_score_treshold] 

    #sample_id = df.sample_id.unique()[0]
    #experiment_id = df.experiment_id.unique()[0]     

    def top3(df):
        df = (df.groupby('ProteinName')['Intensity'].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan)
                  .reset_index())
        return df
    
    dfs = []

    for sample_id in df.sample_id.unique():
        df_sample = df[df.sample_id == sample_id]
        for experiment_id in df_sample.experiment_id.unique():
            
            print(sample_id + "+" + experiment_id)
            df_sample_experiment = df_sample[df_sample.experiment_id == experiment_id]
            df_reduced = df_sample_experiment[["ProteinName", "Intensity"]]
            df_protein = top3(df_reduced)
            #df = df_protein
            specie_mapper = lambda x: x.split("_")[-1]
        
            df_protein["specie"] = df_protein.ProteinName.map(specie_mapper)
            midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
            df_protein = df_protein.set_index(["specie", "ProteinName"])
            res = pd.DataFrame(df_protein.values, columns = midx, index = df_protein.index)
            dfs.append(res)
    df_concat = pd.concat(dfs, axis = 1)
    return df_concat


def top3(df):
    #A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "Pre"].isna().sum(axis=1)<2]
    A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "Pre"].isna().sum(axis=1)<10]
    A = A.iloc[:,A.columns.get_level_values("sample_id") == "Pre"]
    #B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "Post"].isna().sum(axis=1)<2]
    B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "Post"].isna().sum(axis=1)<10]
    B = B.iloc[:,B.columns.get_level_values("sample_id") == "Post"]
    
    print("Finding all overlapping proteins in Sample A and Sample B.")
    # Find overlapping proteins
    overlapping_proteins = list(set(A.index) & set(B.index))
    A = A[A.index.isin(overlapping_proteins)]
    B = B[B.index.isin(overlapping_proteins)]
    
    print("Computing p-values using scipy.stats.ttest_ind() between sample A and B.")
    p_vals = stats.ttest_ind(A, B, axis = 1)[1]
    p_vals = pd.DataFrame(p_vals, columns = ["p"])
    print("Computing q-values.")
    p_vals["q"] = qvalues(p_vals)
    p_vals = pd.DataFrame(p_vals.values, index = A.index, columns = ["p", "q"])
    p_vals.sort_values("q",inplace =True)
    p_vals = p_vals.astype(float)
    
    print("Applying log2-transformation to sample sums")
    A = np.log2(A.sum(axis=1))
    B = np.log2(B.sum(axis=1))
    
    A.name = "1"
    B.name = "2"
    
    df_final = pd.concat([A, B, p_vals], axis = 1)
    print("Computing fold-change")
    df_final["log2(A,B)"] = df_final["1"] - df_final["2"]
    df_final["log2(A,B)"] += 1
    return df_final


df_final["abs(log2)"] = abs(df_final["log2(A,B)"])
result = df_final[(df_final.q < 0.05) & (df_final["abs(log2)"] > 0.58)]
os.chdir("/hdd_14T/data/pxd021197")
result.reset_index().to_csv("top3_ttest.csv", sep = "\t", index = False)

result = pd.read_csv("top3_ttest.csv", sep = "\t")
# 13 differentially abundant proteins.

import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots()
sns.scatterplot(x = df_final["log2(A,B)"], y = -np.log(df_final["q"]), ax = ax)
sns.scatterplot(x = result["log2(A,B)"], y = -np.log(result["q"]), ax = ax)
ax.set_ylabel("-log(q)")
ax.axhline(y = -np.log(0.05), linestyle = "--")
ax.axvline(x = 0.58, linestyle = "--")
ax.axvline(x = -0.58, linestyle = "--")