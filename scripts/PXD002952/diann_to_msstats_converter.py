#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 10:27:31 2021

@author: ptruong
"""

import os 

import pandas as pd
import numpy as np

os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811")

df = pd.read_csv("report.tsv", sep = "\t")
df = df[df["Q.Value"] < 0.01] # filter PSMs Q.Value 166747

df = pd.read_csv("report_recomputed_fdr.tsv", sep = "\t")
df = df[df["fdr"] < 0.01] # filter PSMs Q.Value - 167463



def convert_diann_to_msconvert_aggregated(df):
    df_ = pd.DataFrame()    
    df_["ProteinName"] = df["Protein.Ids"]
    df_["PeptideSequence"] = df["Modified.Sequence"]
    df_["PrecursorCharge"] = df["Precursor.Charge"]
    df_["FragmentIon"] = df["Fragment.Info"]
    df_["ProductCharge"] = np.nan
    df_["IsotopeLabelType"] = "light"
    df_["Intensity"] = df["Fragment.Quant.Raw"] 
    
    bioRep_mapper = lambda x: x.split("Repl")[-1]
    condition_mapper = lambda x: x.split("_")[8]
    run_mapper = lambda x: x.split("_")[5]
    
    df_["BioReplicate"] = df.Run.map(bioRep_mapper)
    df_["Condition"] = df.Run.map(condition_mapper)
    df_["Run"] = df.Run.map(run_mapper)
    return df_

df = convert_diann_to_msconvert_aggregated(df)

# Count and filter based on fragment ion 5-15


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

def filter_n_fragments(df, min_fragments = 4, max_fragments = 6, aggr_fragment_col = "FragmentIon"):
    df["n_fragments"] = df[aggr_fragment_col].str.split(";").apply(len)-1
    df = df[df.n_fragments <= max_fragments]
    df = df[df.n_fragments >= min_fragments]
    df = df.drop("n_fragments", axis = 1)
    return df

df = filter_n_fragments(df, min_fragments = 6, max_fragments = 6, aggr_fragment_col = "FragmentIon")
df = disaggregate(df, fragment_info_col = "FragmentIon", fragment_quant_col = "Intensity", seperator = ";")

df.to_csv("diann_msstats_input_recomputed_fdr_20211201.csv", sep = ",", index = False)

df[df["FragmentIon"] == "b5-unknown^1/398.2033997"]






