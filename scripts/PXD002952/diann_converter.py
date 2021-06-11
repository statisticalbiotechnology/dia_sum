#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 23:08:38 2021

@author: ptruong
"""


import numpy as np
import pandas as pd


def diann_to_triqler(filename):
    df = pd.read_csv(filename, sep = "\t")
    
    run_mapper = lambda x : x.split("_")[5]
    condition_mapper = lambda x : x.split("_")[8]
    
    df["run"] = df["Run"].map(run_mapper)
    df["condition"] = df["Run"].map(condition_mapper)
    df["charge"] = df["Precursor.Charge"]
    #df["searchScore"] = df["CScore"]
    df["searchScore"] = df["Q.Value"]
    df["intensity"] = df['Precursor.Quantity']
    df["peptide"] = df["Stripped.Sequence"]
    df["proteins"] = df["Protein.Ids"]
    
    df_triq = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    return df_triq

df_triq = diann_to_triqler("diann.tsv")
df_triq.to_csv("triqler_input_diann_searchScore_QValue.csv", sep = "\t", index = False)


def diann_to_msstats(filename):
    df = pd.read_csv(filename, sep = "\t")
    run_mapper = lambda x : x.split("_")[5]
    condition_mapper = lambda x : x.split("_")[8]
    replicate_mapper = lambda x : x.split("_")[-1].split("Repl")[1]
    def decoy_mapper(protein):
        if protein.split("_")[0] == "DECOY":
            return True
        else:
            return False
    df['ProteinName'] = df["Protein.Ids"]
    df["decoy"] = df["Protein.Ids"].map(decoy_mapper)
    df = df[df["decoy"] == False]
    df['PeptideSequence'] = df["Stripped.Sequence"]
    df['PrecursorCharge'] = df["Precursor.Charge"]
    df['FragmentIon'] = np.nan
    df['ProductCharge'] = np.nan
    df['IsotopeLabelType'] = "light"
    df['Intensity'] = df["Precursor.Quantity"]
    df["Condition"] = df["Run"].map(condition_mapper)
    df['BioReplicate'] = df["Run"].map(replicate_mapper)
    df["Run"] = df["Run"].map(run_mapper)
    
    df_msstats = df[['ProteinName', 'PeptideSequence', 'PrecursorCharge', 'FragmentIon',
           'ProductCharge', 'IsotopeLabelType', 'Intensity', 'BioReplicate',
           'Condition', 'Run']]
    return df_msstats

df_msstats = diann_to_msstats("diaNN.tsv")
df_msstats.to_csv("msstats_input_diann.tsv", sep = ",", index = False)
