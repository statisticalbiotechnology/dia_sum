#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 08:16:17 2022

@author: ptruong
"""
import os

os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")

import pandas as pd 
import numpy as np

#df = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_002-Pedro_-_Sample_1_-_SW32_-_Repl1.mzML_with_dscore.csv", sep = "\t")
df = pd.read_csv("concatenated_osw_results.csv", sep = "\t")

#m_score filter - m_score at 1% FDR is calculated in R
m_score_treshold = 0.00079433
#df # 5463658
df = df[df.m_score < m_score_treshold] # 82190

def map_zero_transitions(value, threshold):
    """
    Use aggr_Peak_Area to map zero/low intensity transition. 
    """
    x = value.split(";")
    x = np.array(x).astype(float)
    index = np.where(x < threshold)[0]    
    mask = np.ones(x.shape[0], dtype=bool)
    mask[index] = False    
    return mask

def map_top_n_transitions(x, n):    
    x = x.split(";")
    x = np.array(x).astype(float)
    top_n_indices = x.argsort()[-n:][::-1]
    mask = np.zeros(x.shape[0], dtype=bool)
    mask[top_n_indices] = True 
    return mask

def mask_transitions(x, mask):
    x = np.array(x.split(";"))
    return x[mask]


def filter_transition(df, filter_function, *args):
    """
    Function to keep only transition over treshold intensity.
    """
    df["transition_mask"] = df["aggr_Peak_Area"].map(lambda x: filter_function(x, args[0]))
    
    df["aggr_Peak_Area"] = df.apply(lambda x: mask_transitions(x["aggr_Peak_Area"], x["transition_mask"]), axis = 1)
    df["aggr_Peak_Apex"] = df.apply(lambda x: mask_transitions(x["aggr_Peak_Apex"], x["transition_mask"]), axis = 1)
    df["aggr_Fragment_Annotation"] = df.apply(lambda x: mask_transitions(x["aggr_Fragment_Annotation"], x["transition_mask"]), axis = 1)
    df.drop("transition_mask", axis = 1, inplace = True)
    
    df["aggr_Peak_Area"] = df["aggr_Peak_Area"].map(lambda x: ";".join(x))
    df["aggr_Peak_Apex"] = df["aggr_Peak_Apex"].map(lambda x: ";".join(x))
    df["aggr_Fragment_Annotation"] = df["aggr_Fragment_Annotation"].map(lambda x: ";".join(x))

    return df

def convert_osw_to_msstats_aggregated(df):
    df_ = pd.DataFrame()    
    df_["ProteinName"] = df["ProteinName"]
    df_["PeptideSequence"] = df["FullPeptideName"]
    df_["PrecursorCharge"] = df["Charge"]
    df_["FragmentIon"] = df["aggr_Fragment_Annotation"]
    df_["ProductCharge"] = np.nan
    df_["IsotopeLabelType"] = "light"
    df_["Intensity"] = df["aggr_Peak_Area"] 
    
    bioRep_mapper = lambda x: x.split("Repl")[-1].split(".")[0]
    condition_mapper = lambda x: x.split("_")[8]
    run_mapper = lambda x: x.split("_")[5]
    
    df_["BioReplicate"] = df.filename.map(bioRep_mapper)
    df_["Condition"] = df.filename.map(condition_mapper)
    df_["Run"] = df.filename.map(run_mapper)
    return df_


def filter_on_min_peptide(df, n_peptides):
    peptide_count = df.groupby("PeptideSequence").count().ProteinName
    min_peptide = peptide_count[peptide_count >= n_peptides] # greater than
    df = df[df.PeptideSequence.isin(min_peptide.index)]
    return  df

def filter_on_max_peptide(df, n_peptides):
    peptide_count = df.groupby("PeptideSequence").count().ProteinName
    max_peptide = peptide_count[peptide_count <= 10] # less than
    df = df[df.PeptideSequence.isin(max_peptide.index)]
    return df

def drop_decoy_proteins(df):
    df["Decoy"] = (df.ProteinName.map(lambda x:x.split("_")[0]) == "DECOY")
    df = df[df.Decoy != True]
    df.drop("Decoy", inplace = True, axis = 1)
    return df


def disaggregate(df, fragment_info_col = "FragmentIon", fragment_quant_col = "Intensity", seperator = ";"):
    col_fragment_info = fragment_info_col
    col_fragment_quant = fragment_quant_col
    
    fragment_info = df[col_fragment_info].str.split(seperator).apply(pd.Series, 1).stack()
    fragment_info.name = col_fragment_info
    fragment_info.index = fragment_info.index.droplevel(-1)
    fragment_info = fragment_info[fragment_info != ""]
    
    df = df.drop(col_fragment_info, axis = 1)
    
    fragment_quant = df[col_fragment_quant].str.split(seperator).apply(pd.Series, 1).stack()
    fragment_quant.name = col_fragment_quant
    fragment_quant.index = fragment_quant.index.droplevel(-1)
    fragment_quant = fragment_quant[fragment_quant != ""]
    
    df = df.drop(col_fragment_quant, axis = 1)
    
    fragment_quant = pd.DataFrame(fragment_quant)
    fragment_info = pd.DataFrame(fragment_info)
    
    fragment_quant = fragment_quant.reset_index()
    fragment_info = fragment_info.reset_index().drop("index", axis = 1)
    
    fragment_quant[fragment_info_col] = fragment_info
    fragment_quant.set_index("index", inplace=True)
    
    df = df.reset_index()
    df.set_index("index", inplace=True)
    
    df=df.join(fragment_quant)
    df[fragment_quant_col] = df[fragment_quant_col].astype(float)
    return df


#df_ = filter_transition(df, map_zero_transitions, 0.05)  # removes low intensity transitions.
df = filter_transition(df, map_top_n_transitions, 6) # selects top 6 transitions. (Similar to DIA-NN)
df.drop("Unnamed: 0", axis = 1, inplace = True)
df = convert_osw_to_msstats_aggregated(df)
#df # 82190
df = filter_on_min_peptide(df, n_peptides = 2) # 79961
df = filter_on_max_peptide(df, n_peptides = 10) # 62539
df = filter_n_fragments(df, min_fragments = 1, max_fragments = 6, aggr_fragment_col = "FragmentIon")
df = drop_decoy_proteins(df) # 62350
´
df = disaggregate(df, fragment_info_col = "FragmentIon", fragment_quant_col = "Intensity", seperator = ";")

df.to_csv("osw_msstats_input.csv", sep = ",", index = False)


#df_.to_csv("concatenated_osw_results_transitions_filtered_n_6.csv", sep = "\t", index = False)































