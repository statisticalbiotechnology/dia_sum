#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:54:17 2022

@author: ptruong
"""

import pandas as pd
import numpy as np

def parse_triqler(triqler_output_file):
    """
    Parses triqler output format to pandas dataframe.
    """
    f = open(triqler_output_file, "r")
    lines = f.readlines()
    line = lines.pop(0)
    cols = line.split("\n")[0].split("\t")[:]
    n_cols = len(cols)
    
    data_array = []
    for line in lines:
        line = line.split("\n")[0].split("\t")
        vals = line[:n_cols-1]
        peptides = ";".join(line[n_cols-1:])
        data = vals + [peptides]
        data_array.append(data)
    df = pd.DataFrame(data_array, columns = cols)

    df = pd.concat([df[["protein", "peptides"]], df.drop(["protein", "peptides"], axis = 1).astype(float)], axis = 1)
    
    return df



lib_pxd = pd.read_csv("lib.tsv", sep = "\t")
lib_input = pd.read_csv("library_specie.tsv", sep = "\t")
            
lib_pxd = lib_pxd[["PrecursorMz", "ProductMz", "ProteinName", "Genes", "PeptideSequence", "ModifiedPeptide",
         "PrecursorCharge", "LibraryIntensity", "Tr_recalibrated", "UniprotID", "FragmentType",
         "FragmentCharge", "FragmentSeriesNumber", "FragmentLossType"]]

lib_pxd.rename({"ProteinName":"ProteinId", "Genes":"GeneName", "ModifiedPeptide":"ModifiedPeptideSequence",
                "Tr_recalibrated":"NormalizedRetentionTime"}, axis = 1, inplace = True)

lib_pxd["PrecursorIonMobility"] = np.nan

lib_pxd["Annotation"] = np.nan

lib_pxd = lib_pxd[["PrecursorMz", "ProductMz", "Annotation", "ProteinId", "GeneName", "PeptideSequence", 
         "ModifiedPeptideSequence", "PrecursorCharge", "PrecursorIonMobility", "LibraryIntensity", 
         "NormalizedRetentionTime", "UniprotID", "FragmentType", "FragmentCharge", "FragmentSeriesNumber",
         "FragmentLossType"]]




lib_pxd.to_csv("lib_converted.tsv", sep = "\t", index = False)

#### Fix - on added decoy

import os

os.chdir("/hdd_14T/data/PXD029721_PBMC_dia_sum_test/ftp.pride.ebi.ac.uk/pride/data/archive/2022/04/PXD029721/test/lib_converted_decoy")

lib = pd.read_csv("lib_converted_decoy.tsv", sep  = "\t")
lib = pd.read_csv("2021-9_fullExp-lib.decoy.tsv", sep = "\t")
lib = lib[["PrecursorMz", "ProductMz", "ProteinId", "GeneName", "PeptideSequence", "ModifiedPeptideSequence",
     "PrecursorCharge", "LibraryIntensity", "NormalizedRetentionTime",
     "PrecursorIonMobility", "FragmentType", "ProductCharge", 
     "FragmentSeriesNumber"]]

lib.rename({"ProductCharge": "FragmentCharge"}, axis = 1, inplace = True)
lib.to_csv("library_converted_decoy_reformatted.tsv", sep = "\t", index = False)


## ABVOE AWORKS 
lib[lib.PeptideSequence == "GVEKPSDSLLLQAPDGYR"]
lib[lib.PeptideSequence == "LLADELISPK"]
d = lib[lib.ProteinId == "DECOY_KPBB_HUMAN"].PeptideSequence.unique()
t = lib[lib.ProteinId == "KPBB_HUMAN"].PeptideSequence.unique()
list(set(t) & set(d))
lib[lib.GeneName == "PHKB"].ProteinId.unique()
lib[lib.GeneName == "PHKB"].PeptideSequence.unique()
 

### Check decoy on LFQBENCH

os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707")

lib_lfq = pd.read_csv("library_specie_decoy_formatted.tsv", sep = "\t")

# checl TKHVTIEIYNLDTEQSR
lib_lfq[lib_lfq.PeptideSequence == "TKHVTIEIYNLDTEQSR"].ProteinId
t = lib_lfq[lib_lfq.ProteinId == "Q9NSE4_HUMAN"].PeptideSequence.unique()
d = lib_lfq[lib_lfq.ProteinId == "DECOY_Q9NSE4_HUMAN"].PeptideSequence.unique()
list(set(t) & set(d))
lib_lfq[lib_lfq.GeneName == "IARS2"].ProteinId.unique()
lib_lfq[lib_lfq.GeneName == "IARS2"].PeptideSequence.unique()

"""
df = pd.read_csv("PBMC_5000_DIA-NN_GPF_based_report.tsv", sep = "\t")
len(df["Protein.Group"].unique())


len(df[df["Protein.Q.Value"] < 0.01]["Protein.Group"].unique())
len(df[df["Protein.Q.Value"] < 0.01]["Protein.Group"].unique())
len(df[df["Q.Value"] < 0.01]["Protein.Group"].unique())

len(df[df["Protein.Q.Value"] < 0.05]["Protein.Group"].unique())
len(df[df["Protein.Q.Value"] < 0.05]["Protein.Group"].unique())
len(df[df["Q.Value"] < 0.05]["Protein.Group"].unique())


len(df[(df["Q.Value"] < 0.01) & (df["PG.Q.Value"] < 0.01)]["Protein.Group"].unique())
len(df[(df["Q.Value"] < 0.05) & (df["PG.Q.Value"] < 0.05)]["Protein.Group"].unique())

len(df[(df["Q.Value"] < 0.05) & (df["PG.Q.Value"] < 0.05)]["Protein.Names"].unique())



df.columns
"""
triqler results
triqler = parse_triqler("proteins.tsv")

