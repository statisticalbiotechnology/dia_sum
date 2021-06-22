#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 12:08:41 2021

@author: ptruong
"""

import numpy as np
import pandas as pd 

import os 
# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_to_msqrobsum_converter import get_mSqRobSum_input
from get_protein_specie_map_from_fasta import fasta_to_protein_specie_map
# db 
os.chdir("/home/ptruong/git/dia_sum/database/2021-06-16")

fasta = "no_shared_HUMAN_ECOLI_YEAST.fasta"
protein_specie_map = fasta_to_protein_specie_map(fasta)

# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire/msfragger")

# Add specie to Easypqp format.
def append_specie_to_easyPQP_library(filename = "library.tsv"):
    lib = pd.read_csv(filename, sep = "\t")
    lib["specie"] = lib.ProteinId.map(protein_specie_map)
    lib["ProteinId"] = lib["ProteinId"] + "_" + lib["specie"]
    lib.drop(["specie"], axis = 1, inplace = True)
    lib.to_csv(f"{filename[:-4]}_specie.tsv", sep = "\t", index = False)

append_specie_to_easyPQP_library(filename = "library.tsv")

# Run OpenSwathDecoyGenerator

# Convert OpenSwathDecoyGenerator output to EasyPQP format.
def format_openSwathDecoyGenerator_to_diann(filename = "library_specie_decoy.tsv"):
    lib_decoy = pd.read_csv(filename, sep = "\t")
    lib_decoy["FragmentCharge"] = lib_decoy["ProductCharge"]
    #lib_decoy = pd.concat([lib_decoy.iloc[:, lib_decoy.columns.isin(lib.columns)], lib_decoy["FragmentCharge"]], axis = 1)
    lib_decoy["FragmentLossType"] = np.nan
    lib_decoy = lib_decoy[lib.columns.unique()]
    # Add as many column as possible 
    lib_decoy.to_csv(f"{filename[:-4]}_formatted.tsv", sep = "\t", index=False)

format_openSwathDecoyGenerator_to_diann(filename = "library_specie_decoy.tsv")

