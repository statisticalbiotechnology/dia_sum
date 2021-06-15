#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:14:25 2021

@author: ptruong
"""


import os 
# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from diann_to_triqler_converter import diann_to_triqler
from triqler_to_msqrobsum_converter import get_mSqRobSum_input
from get_protein_specie_map_from_fasta import fasta_to_protein_specie_map
# db 
os.chdir("/home/ptruong/git/dia_sum/database")
protein_specie_map = fasta_to_protein_specie_map("napedro_3mixed_human_yeast_ecoli_20140403_iRT_reverse.fasta")

# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/res_20210609_DIAUmpire/MSFragger/diann_reformatted_lib")

# Convert diann_to_triqler

df_triq = diann_to_triqler("diann.tsv")
df_triq.to_csv("triqler_input_diann.csv", sep = "\t", index = False)

triqler = pd.read_csv("triqler_input_diann.csv", sep = "\t")

triqler.dropna(inplace = True) # drop rows with "nan" protein name
# Fix protein labelling

def remove_decoy_tag(protein):
    if protein.split("_")[0] == "DECOY":
        return protein.split("_")[1]
    else:
        return protein.split("_")[0]
    
triqler["protein_nonTagged"] = triqler.proteins.map(remove_decoy_tag)
triqler["specie"] = triqler["protein_nonTagged"].map(protein_specie_map)
triqler["proteins"] = triqler["proteins"] + "_" + triqler["specie"]
triqler.drop(["protein_nonTagged", "specie"], axis = 1, inplace=True)
expr_, fd_, pd_ = get_mSqRobSum_input(triqler)

os.chdir("/hdd_14T/data/PXD002952/res_20210604_DIAUmpire/MSFragger/diann/msqrob")
remap_sample = lambda s : "X"+".".join(s.split("-"))
pd_["sample"] = pd_["sample"].map(remap_sample)

expr_.to_csv("expr.csv", sep = "\t", index = False)
fd_.to_csv("fd.csv", sep = "\t", index = False)
pd_.to_csv("pd.csv", sep = "\t", index = False)



# Relabel library protein
library = pd.read_csv("library.tsv")
