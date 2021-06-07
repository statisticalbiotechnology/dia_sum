#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 18:06:32 2021

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
os.chdir("/home/ptruong/git/dia_sum/database")
protein_specie_map = fasta_to_protein_specie_map("napedro_3mixed_human_yeast_ecoli_20140403_iRT_reverse.fasta")

# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/res_20210604_DIAUmpire/MSFragger")

lib = pd.read_csv("library.tsv", sep = "\t")
lib["specie"] = lib.ProteinId.map(protein_specie_map)
lib["ProteinId"] = lib["ProteinId"] + "_" + lib["specie"]
lib.drop(["specie"], axis = 1, inplace = True)
lib.to_csv("library_specie.tsv", sep = "\t", index = False)