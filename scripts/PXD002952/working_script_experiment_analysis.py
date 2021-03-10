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

df_t = df_merged[df_merged.m_score < 0.01] #FDR-treshold
df = df_t[df_t.decoy == 0]

#top 3 peptides to select protein

def get_protein_quantity_df(df_sample):
    """
    E.g. input of df_sample:
        
    df_sample = df[df.experiment_id == "001-Pedro"]
    """
    def get_protein_quantity(protein):
        df_prot = df_sample[df_sample.ProteinName == protein]
        if len(df_prot) > 2:
            intensity = df_prot.Intensity.sort_values(ascending=False).head(3).mean()
        else:
            intensity = np.nan
        return intensity
    
    import time
    start=time.time()
    df_sample["protein_quantity"] = df_sample["ProteinName"].map(get_protein_quantity)
    end=time.time()
    print("time: " + str(end-start) + " s")
    prot = df_sample[df_sample.protein_quantity == 199.523][["experiment_id", "sample_id", "ProteinName","protein_quantity"]]
    prot.drop_duplicates()
    prot = df_sample[df_sample["protein_quantity"].notna()]
    prot = prot[["experiment_id", "sample_id", "ProteinName", "protein_quantity"]]
    prot = prot.drop_duplicates()
    return prot


def protein_quantity_df(df):
    """
    df - df_merged, tresholded and removed decoy protein
    
    
    """
    df_protein = pd.DataFrame()
    for sample_id in df.experiment_id.unique():
        df_sample = df[df.experiment_id == sample_id]    
        df_protein_subset = get_protein_quantity_df(df_sample)
        df_protein = df_protein.append(df_protein_subset)
    return df_protein
# Summary

exp_array = []
sample_array = []
proteins_array = []
for exp in df_protein.experiment_id.unique():
    df_exp = df_protein[df_protein.experiment_id == exp]
    sample = df_exp.sample_id.unique()
    n_proteins = len(df_exp)
    exp_array.append(exp)
    sample_array.append(sample)
    proteins_array.append(n_proteins)

df_summary = pd.DataFrame(np.array([exp_array, sample_array, proteins_array]).T, columns = ["experiment_id", "sample_id", "proteins"])

























