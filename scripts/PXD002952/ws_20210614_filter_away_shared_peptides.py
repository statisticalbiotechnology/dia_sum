#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 01:39:05 2021

@author: ptruong


It was messy the randomly select proteins with shared peptide to retain some proteins because
some shared peptides are among many proteins. I just filter away all the shared proteins.

"""


import os 
import re
import pandas as pd
import numpy as np

# read in .fasta file and count shared peptides

os.chdir("/home/ptruong/git/dia_sum")
filename = "database/2021-06-07/UP00000625_UP000002311_UP000005640_reformatted.fasta"
#filename = "database/napedro_3mixed_human_yeast_ecoli_20140403_iRT_reverse.fasta"

file = open(filename, "r")

protein_list = []
sequence_list = []
for line in file: 
    if line[0] == ">":
        protein = line 
    else:
        sequence = line.rstrip()
        split_sequence = re.split(r"(?<=[KR])", sequence)
        split_sequence = list(dict.fromkeys(split_sequence))
        sequence_list += split_sequence
        protein_list += [protein for i in range(len(split_sequence))]
        

df_ = pd.DataFrame(np.array([protein_list, sequence_list]).T, columns = ["protein", "sequence"])
df = pd.DataFrame(np.array([protein_list, sequence_list]).T, columns = ["protein", "sequence"])

def decoy_map(protein):
    if protein.split("_")[0] == ">reverse":
        return True
    else:
        return False
    
    
df["decoy"] = df.protein.map(decoy_map)


df = df[df.decoy == False]
df["seq_length"] = df.sequence.str.len()
df = df[df["seq_length"] > 7]
df.drop("seq_length", axis = 1, inplace = True)
df.drop_duplicates(inplace=True)
counted_df = df.groupby("sequence").count().sort_values(by = "protein", ascending = False)

# We must reverse the thinking - find protein names to filter away from original .fasta

multi_protein_seq = counted_df[counted_df > 1].dropna().index # sequences that belong to multiple proteins.

len(df[df.sequence.isin(multi_protein_seq)].protein.unique())

len(df_.protein.unique())

filter_list = df[df.sequence.isin(multi_protein_seq)].protein.unique() #proteins to filter away from original .fasta

file = open("database/UP00000625_UP000002311_UP000005640_reformatted.fasta", "w")
for i in range(len(df_)):
    elem = df_.iloc[i]
    protein = elem.protein
    if protein not in filter_list:
        file.write(elem.protein)
        file.write(elem.sequence + "\n")
    

    
    
    
    
    
    
    
    
    
    