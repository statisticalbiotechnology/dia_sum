#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 15:46:08 2021

@author: ptruong
"""

import os 

import pandas as pd
import numpy as np

# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from get_protein_specie_map_from_fasta import fasta_to_protein_specie_map
#db
os.chdir("/home/ptruong/git/dia_sum/database")
protein_specie_map = fasta_to_protein_specie_map("2021-04-27-decoys-reviewed-contam-UP000005640-UP000002311-UP000000625.fas")

#### TEST with OSW ######
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")
ms_osw = pd.read_csv("msstats_20210511_unfiltered.csv", sep = ",")



#### START HERE ##### 20210603
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")
ms_diann = pd.read_csv("msstats_input_diann.tsv", sep = ",")



len(ms_diann.ProteinName.unique()) - len(ms_osw.ProteinName.unique())

def decoy_map(protein):
    if protein.split("_")[0] == "DECOY":
        return True
    else:
        return False
    
    
ms_diann["decoy"] = ms_diann.ProteinName.map(decoy_map)
ms_osw["decoy"] = ms_osw.ProteinName.map(decoy_map)

ms_osw.decoy.sum()
ms_diann.decoy.sum() # Why is there no decoy here?



# diann --> aligned format --> msstat via swath2stat converter <-------------
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")
aligned = pd.read_csv("aligned.csv", sep = "\t")

os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")
df = pd.read_csv("diaNN.tsv", sep = "\t")



ProteinName 
PeptideSequence 
PrecursorCharge 
Condition 
BioReplicate 
Run      
RT 
NakedSequence 
FragmentIon 
Intensity

def remove_decoy_tag(protein):
    if protein.split("_")[0] == "DECOY": 
        return protein.split("_")[1]
    else:
        return protein.split("_")[0]

def decoy_map(protein):
    if protein.split("_")[0] == "DECOY":
        return True
    else:
        return False
    
df["protein_nonTagged"] = df["Protein.Ids"].map(remove_decoy_tag)
df["specie"] = df["protein_nonTagged"].map(protein_specie_map)
df["Protein.Ids"] = df["Protein.Ids"] + "_" + df["specie"]
df["decoy"] = df["Protein.Ids"].map(decoy_map)
df["ProteinName"] = df["Protein.Ids"]

df["PeptideSequence"] = df["Precursor.Id"]
df["PrecursorCharge"] = df["Precursor.Charge"]


condition_mapper = lambda x : x.split("_")[8]
replicate_mapper = lambda x : x.split("_")[-1].split("Repl")[-1]
run_mapper = lambda x: x.split("_")[5]

df["Condition"] = df.Run.map(condition_mapper)
df["BioReplicate"] = df.Run.map(replicate_mapper)
df["Run"] = df.Run.map(run_mapper)
df["RT"] = df["RT"]
df["NakedSequence"] = df["Stripped.Sequence"]
df["FragmentIon"] = 
df["Intensity"] = 


len(aligned.Sequence.unique()) # 508917 


aligned




## 

os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")

msqrob = pd.read_csv("msqrobsum_result.csv", sep = "\t")




