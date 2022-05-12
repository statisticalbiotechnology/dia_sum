#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 19:06:44 2022

@author: ptruong
"""
import os
import pandas as pd 
import numpy as np
import argparse


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


def filter_and_convert_osw_to_msstats(input_file, output, m_score_threshold_file, n_transitions = 6, min_peptides = 2, max_peptides = 10):
    m_score_threshold = float(open(m_score_threshold_file).read())
    df = pd.read_csv(input_file, sep = "\t")
    df = df[df.m_score < m_score_threshold] # 82190
    #print(df)
    df = filter_transition(df, map_top_n_transitions, n_transitions) # selects top 6 transitions. (Similar to DIA-NN)
    df.drop("Unnamed: 0", axis = 1, inplace = True)
    df = convert_osw_to_msstats_aggregated(df)
    df = filter_on_min_peptide(df, n_peptides = min_peptides) 
    df = filter_on_max_peptide(df, n_peptides = max_peptides) 
    df = drop_decoy_proteins(df) 
    
    df = disaggregate(df, fragment_info_col = "FragmentIon", fragment_quant_col = "Intensity", seperator = ";")
    df.drop_duplicates(inplace = True)
    df = df.groupby(["PeptideSequence", "FragmentIon", "Run"]).max()
    df = df.reset_index()
    df.to_csv(output, sep = ",", index = False)

    

parser = argparse.ArgumentParser(
    description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots calibration plots.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='Input osw file. This should be the concatenated file.')

parser.add_argument('--output', type=str,
                    help='Output name.')

parser.add_argument('--m_score_threshold_file', type=str,
                    help='mscore threshold. It should be computed with mscore4pepfdr.')

parser.add_argument('--n_transitions', type=int, default = 6,
                    help='Maximum transitions.')

parser.add_argument('--min_peptides', type=int, default = 2,
                    help='Min number of peptides to allow for a protein.')

parser.add_argument('--max_peptides', type=int, default = 10,
                    help='Max number of peptides to allow for a protein.')


# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
output = args.output
m_score_threshold_file = args.m_score_threshold_file
n_transitions = args.n_transitions
min_peptides = args.min_peptides
max_peptides = args.max_peptides
                      

#m_score_treshold = 0.00079433
#n_transitions = 6
#min_peptides = 2
#max_peptides = 10
#input_file = "concatenated_osw_results.csv"
#output = "osw_msstats_input.csv"

if __name__ == "__main__":
    print(output)
    filter_and_convert_osw_to_msstats(input_file = input_file, output = output, m_score_threshold_file = m_score_threshold_file, n_transitions = n_transitions, min_peptides = min_peptides, max_peptides = max_peptides)



