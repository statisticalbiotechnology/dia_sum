#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 20:56:18 2021

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

# Filtering
aligned = aligned[aligned.m_score < 0.01] # q-value treshold
aligned = aligned[aligned.decoy != True] # remove decoys


def get_protein_intensity_df(aligned, run):
    """
    aligned - aligned.df
    run - e.g. - run_list[0]
    """
    df_run = aligned[aligned["run"] == run]
    protein_list = df_run.ProteinName.unique()
    
    def get_protein_intensity(df_run, protein):
        df_protein = df_run[df_run.ProteinName == protein]
        if len(df_protein) > 2:
            intensity = df_protein.Intensity.sort_values(ascending=False).head(3).mean()
        else:
            intensity = np.nan
        return intensity
    
    protein_array = []
    for protein in protein_list:
        intensity = get_protein_intensity(df_run, protein)
        protein_array.append([protein, intensity])
        
    protein_intensity_df = pd.DataFrame(protein_array, columns = ["protein", run]).set_index("protein")
    
    sample = "_".join(list(df_run["sample"].unique())) #using join to see if this messes up...
    res = {sample: protein_intensity_df}
    
    #protein_intensity_df["run"] = run
    return res

def get_all_sample_protein_intensities_no_sample_info(aligned):
    """
    aligned - filtered aligned dataframe.
    """
    df_res = pd.DataFrame()
    for run in run_list:
        protein_intensity_df = get_protein_intensity_df(aligned, run)
        df_res = pd.concat([df_res, protein_intensity_df], axis = 1)
    return df_res

def get_all_sample_protein_intensities(aligned):
    """
    aligned - filtered aligned dataframe.
    
    Hierarchical dataframe with sample information.
    
    Top level - sample
    Low level - run
    
    """
    res = {}
    for run in run_list:
        print(run)
        partial_res = get_protein_intensity_df(aligned, run)
        if list(partial_res.keys())[0] in list(res.keys()):
            sample = list(partial_res.keys())[0]
            res = {**res, **{sample : pd.concat([res[sample], partial_res[sample]], axis=1)}}
        else:
            res = {**res, **partial_res}
        
    df_res = pd.concat(res, axis = 1)
    df_res["species"] = df_res.index.map(species_map) # Add species col
    return df_res

###

df = get_all_sample_protein_intensities(aligned)







###########
# Triqler #
###########

from triqler_output_to_df import parse_triqler

triqler = parse_triqler("proteins.tsv")

triqler[triqler.q_value < 0.05]
triqler[triqler.posterior_error_prob < 0.05]
