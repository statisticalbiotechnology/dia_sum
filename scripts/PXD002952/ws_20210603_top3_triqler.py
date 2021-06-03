#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:42:19 2021

@author: ptruong
"""

import os 
import pandas as pd
import numpy as np
# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_output_to_df import parse_triqler

# db
from get_protein_specie_map_from_fasta import fasta_to_protein_specie_map
os.chdir("/home/ptruong/git/dia_sum/database")
protein_specie_map = fasta_to_protein_specie_map("2021-04-27-decoys-reviewed-contam-UP000005640-UP000002311-UP000000625.fas")
#protein_specie_map = protein_specie_map.set_index("protein").T.to_dict()
# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")

df = pd.read_csv("diaNN.tsv", sep = "\t")
triq = pd.read_csv("triqler_input_diann.csv", sep = "\t")                      

def compute_triqler_top3_submodule(run):
    triq_run = triq[triq.run == run]                    
    def triqler_top3(triq_run):
        res = triq_run.groupby("proteins")["intensity"].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan).reset_index()
        return res
    
    def triqler_printout_unique_peptides_proteins(run):
        triq_run = triq[triq.run == run]
        condition = triq_run.condition.unique()
        print(f"run : {run} - condition : {condition}")
        print(f"Unique peptides detected: {len(triq_run.peptide.unique())}")
        print(f"Unique proteins detected: {len(triq_run.proteins.unique())}")
        print()
    
    triqler_printout_unique_peptides_proteins((run))
    
    res = triqler_top3(triq_run)
    def remove_decoy_tag_protein(protein):
        if protein.split("_")[0] == "DECOY":
            return protein.split("_")[1]
        else: return protein.split("_")[0]
    
    res["proteins_nonTagged"] = res.proteins.map(remove_decoy_tag_protein)
    res["specie"] = res.proteins_nonTagged.map(protein_specie_map)
    experiment_id = triq_run.run.unique()[0]
    sample_id = triq_run.condition.unique()[0]
    res["ProteinName"] = res["proteins"]
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    res = res.set_index(["specie", "ProteinName"])
    res = res.drop(["proteins", "proteins_nonTagged"], axis = 1)
    df = pd.DataFrame(res.values, columns = midx, index = res.index)
    return df


df = compute_triqler_top3_submodule(run)

dfs = []
for run in triq.run.unique():
    dfs.append(compute_triqler_top3_submodule(run))
     
df = pd.concat(dfs, axis = 1)
df # Check why we have more species than in the fasta??? Is this because of how we build the spectral library?
df = df[df.index.get_level_values("specie").isin(["HUMAN", "ECOLI", "YEAST"])]

df[df.index.get_level_values("specie").isin(["HUMAN"])]
df[df.index.get_level_values("specie").isin(["ECOLI"])]
df[df.index.get_level_values("specie").isin(["YEAST"])]







