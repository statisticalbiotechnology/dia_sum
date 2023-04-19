#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 20:17:16 2021

@author: ptruong
"""

import pandas as pd 
import numpy as np

aligned = pd.read_csv("aligned.csv", sep = "\t")

species_map = lambda x : x.split("_")[-1]
run_map = lambda x : x.split("_")[5]
sample_map = lambda x : x.split("_")[8]


species_list = aligned.ProteinName.map(species_map).unique()
run_list = aligned.filename.map(run_map).unique()
sample_list = aligned.filename.map(sample_map).unique()

aligned["specie"] = aligned.ProteinName.map(species_map)
aligned["run"] = aligned.filename.map(run_map)
aligned["sample"] = aligned.filename.map(sample_map)

# treshold

def multiIndex_aligned(aligned):
    tuples = [list(aligned.index), list(aligned["specie"]), list(aligned["run"]), list(aligned["sample"])]
    tuples = list(zip(*tuples))
    midx = pd.MultiIndex.from_tuples(tuples, names=["id", "specie", "run", "sample"])
    aligned = pd.DataFrame(aligned.drop(["specie", "run", "sample"],axis=1).values, 
                           columns = aligned.drop(["specie", "run", "sample"], axis=1).columns, index=midx)
    return aligned

aligned = multiIndex_aligned(aligned)

aligned.columns

# Filter 
aligned = aligned[aligned["decoy"] == 0]
aligned = aligned[aligned["m_score"] < 0.05]

run = aligned.index.get_level_values("run").unique()[0]
df_run = aligned[aligned.index.get_level_values("run") == run] 

def get_protein_intensity(df_run, protein):
    df_protein = df_run[df_run.ProteinName == protein]
    if len(df_protein) > 2:
        intensity = df_protein.Intensity.sort_values(ascending=False).head(3).mean()
    else:
        intensity = np.nan
    return intensity


protein_list = df_run.ProteinName.unique()
protein_array = []
for protein in protein_list:
    intensity = get_protein_intensity(df_run, protein)
    protein_array.append([protein, intensity])



