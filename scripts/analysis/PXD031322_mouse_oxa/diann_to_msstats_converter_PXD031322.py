#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:28:27 2022

@author: ptruong
"""

import os 

import pandas as pd
import numpy as np
import argparse


def convert_diann_to_msconvert_aggregated(df):
    df_ = pd.DataFrame()    
    df_["ProteinName"] = df["Protein.Ids"]
    df_["PeptideSequence"] = df["Modified.Sequence"]
    df_["PrecursorCharge"] = df["Precursor.Charge"]
    df_["FragmentIon"] = df["Fragment.Info"]
    df_["ProductCharge"] = np.nan
    df_["IsotopeLabelType"] = "light"
    df_["Intensity"] = df["Fragment.Quant.Raw"] 
    
    #bioRep_mapper = lambda x: x.split("Repl")[-1]
    #condition_mapper = lambda x: x.split("_")[8]
    run_mapper = lambda x: x.split("_")[5]
    
    df_["BioReplicate"] = df.Run.map(lambda x:x.split("-")[-1][-1])
    df_["Condition"] = df.Run.map(lambda x:x.split("-")[-1][:-1])
    df_["Run"] = df.Run.map(lambda x:x.split("-")[-1])
    df_=df_[df_.Condition.isin(["Ctrl", "ST", "LT"])] # Keep only Ctrl, ST, LTs
    return df_

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


def filter_on_min_peptide(df, n_peptides):
    peptide_count = df.groupby("ProteinName").count().PeptideSequence
    min_peptide = peptide_count[peptide_count >= n_peptides] # greater than
    df = df[df.ProteinName.isin(min_peptide.index)]
    return  df

def filter_on_max_peptide(df, n_peptides):
    peptide_count = df.groupby("ProteinName").count().PeptideSequence
    max_peptide = peptide_count[peptide_count <= n_peptides] # less than
<<<<<<< Updated upstream
    df = df[df.ProteinName.isin(max_peptide.index)]
=======
    df = df[df.PeptideSequence.isin(max_peptide.index)]
>>>>>>> Stashed changes
    return df

def drop_decoy_proteins(df):
    df["Decoy"] = (df.ProteinName.map(lambda x:x.split("_")[0]) == "DECOY")
    df = df[df.Decoy != True]
    df.drop("Decoy", inplace = True, axis = 1)
    return df



def convert_diann_to_msstats(input_file = "report_recomputed_fdr.tsv", output = "msstats_input.csv", fdr_threshold = 0.01, qCol = "q"):
    print("Converting diann to MSstats")
    print(f"fdr threshold {fdr_threshold}")
    df = pd.read_csv(input_file, sep = "\t")
    df = df[df[qCol] < fdr_threshold] # filter PSMs Q.Value - 167463

    df = convert_diann_to_msconvert_aggregated(df)
        
    #df # 167463
    df = filter_on_min_peptide(df, n_peptides = 2) # 166236
    df = filter_on_max_peptide(df, n_peptides = 14) # 93395
    df = filter_n_fragments(df, min_fragments = 1, max_fragments = 6, aggr_fragment_col = "FragmentIon")
    df = drop_decoy_proteins(df)
    
<<<<<<< Updated upstream
    #count_table = df.groupby("PeptideSequence")["Condition"].unique()
    #count_table = (count_table.apply(len) > 2)
    #count_table = count_table[count_table == True]
    
    #df = df[df.PeptideSequence.isin(count_table.index)]
=======
    count_table = df.groupby("PeptideSequence")["Condition"].unique()
    count_table = (count_table.apply(len) > 2)
    count_table = count_table[count_table == True]
    
    df = df[df.PeptideSequence.isin(count_table.index)]
>>>>>>> Stashed changes
    df = disaggregate(df, fragment_info_col = "FragmentIon", fragment_quant_col = "Intensity", seperator = ";")
    
    df.to_csv(output, sep = ",", index = False)
    df[df.Condition.isin(["Ctrl", "ST"])].to_csv("Ctrl_ST_"+output, sep = ",", index = False)
    df[df.Condition.isin(["Ctrl", "LT"])].to_csv("Ctrl_LT_"+output, sep = ",", index = False)
    df[df.Condition.isin(["LT", "ST"])].to_csv("LT_ST_"+output, sep = ",", index = False)
    df.to_csv("LT_ST_Ctrl_"+output, sep = ",", index = False)
    


df_red = df.copy()
df_full = df.copy()

len(df_red.ProteinName.unique())
len(df_full.ProteinName.unique())

parser = argparse.ArgumentParser(
    description='Convert DIA-NN to MSstats input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--fdr_threshold', type=float, default = 0.01,
                    help='fdr threshold to apply on qCol.')

parser.add_argument('--qCol', type=str, default = "q",
                    help='q value column to use.')

parser.add_argument('--output', type=str, default = "msstats_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
fdr_threshold = args.fdr_threshold
qCol = args.qCol
output = args.output


if __name__ == "__main__":
    convert_diann_to_msstats(input_file, output, fdr_threshold, qCol)



