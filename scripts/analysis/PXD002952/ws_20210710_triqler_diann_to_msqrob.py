#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:14:25 2021

@author: ptruong
"""


import os 
# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_to_msqrobsum_converter import get_mSqRobSum_input
from get_protein_specie_map_from_fasta import fasta_to_protein_specie_map
# db 
os.chdir("/home/ptruong/git/dia_sum/database")
protein_specie_map = fasta_to_protein_specie_map("2021-04-27-decoys-reviewed-contam-UP000005640-UP000002311-UP000000625.fas")

# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")

triqler = pd.read_csv("triqler_input_diann.csv", sep = "\t")
# Fix protein labelling

def remove_decoy_tag(protein):
    try:
        if protein.split("_")[0] == "DECOY":
            return protein.split("_")[1]
        else:
            return protein.split("_")[0]
    except:
        print(protein)
triqler = triqler[~triqler.proteins.isna()]
triqler["protein_nonTagged"] = triqler.proteins.map(remove_decoy_tag)
triqler["specie"] = triqler["protein_nonTagged"].map(protein_specie_map)
triqler["proteins"] = triqler["proteins"] + "_" + triqler["specie"]
triqler.drop(["protein_nonTagged", "specie"], axis = 1, inplace=True)
expr_, fd_, pd_ = get_mSqRobSum_input(triqler)
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger/msqrob_input")
remap_sample = lambda s : "X"+".".join(s.split("-"))
pd_["sample"] = pd_["sample"].map(remap_sample)

expr_.to_csv("expr.csv", sep = "\t", index = False)
fd_.to_csv("fd.csv", sep = "\t", index = False)
pd_.to_csv("pd.csv", sep = "\t", index = False)
