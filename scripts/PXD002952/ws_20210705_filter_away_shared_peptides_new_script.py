#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  5 12:32:50 2021

@author: ptruong
"""


import os 
import re
import pandas as pd
import numpy as np

# read in .fasta file and count shared peptides

os.chdir("/home/ptruong/git/dia_sum")
filename = "database/2021-06-07/UP00000625_UP000002311_UP000005640.fasta"

#filename = "no_shared_HUMAN_ECOLI_YEAST.fasta"
#filename = "database/napedro_3mixed_human_yeast_ecoli_20140403_iRT_reverse.fasta"

file = open(filename, "r")


#line = file.readline()
protein_list = []
sequence_list = []
sequence = ""
for line in file: 
    if line[0] == ">":
        if len(sequence) > 0:
            split_sequence = re.split(r"(?<=[KR])", sequence)
            split_sequence = list(dict.fromkeys(split_sequence))
            sequence_list += split_sequence
            protein_list += [protein for i in range(len(split_sequence))]
        protein = line
        sequence = ""
    else:
        sequence += line.rstrip()
#        split_sequence = re.split(r"(?<=[KR])", sequence)
#        split_sequence = list(dict.fromkeys(split_sequence))
#        sequence_list += split_sequence
#        protein_list += [protein for i in range(len(split_sequence))]
        

#df_ = pd.DataFrame(np.array([protein_list, sequence_list]).T, columns = ["protein", "sequence"])
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

file = open("UP00000625_UP000002311_UP000005640_no_shared_peptides.fasta", "w")

count = 0
read_file = open(filename, "r")

for line in read_file: 
    if line[0] == ">":
        protein = line
        if protein not in filter_list:
            write_sequence = True
            file.write(protein)
        else:
            write_sequence = False
    else:
        if write_sequence == True:
            sequence = line
            file.write(sequence)
        else:
            continue


    
filename = "UP00000625_UP000002311_UP000005640_no_shared_peptides.fasta"
    
# Count proteins

os.chdir("/home/ptruong/git/dia_sum")
filename = "database/2021-06-07/UP00000625_UP000002311_UP000005640.fasta"

def read_fasta_to_df(filename):
    file = open(filename, "r")
    protein_list = []
    sequence_list = []
    sequence = ""
    for line in file: 
        if line[0] == ">":
            if len(sequence) > 0:
                sequence_list.append(sequence)
                protein_list.append(protein)
            protein = line
            sequence = ""
        else:
            sequence += line.rstrip()
    df = pd.DataFrame(np.array([protein_list, sequence_list]).T, columns = ["protein", "sequence"])
    return df    
        
df_raw = read_fasta_to_df(filename)

filename = "database/2021-07-05/UP00000625_UP000002311_UP000005640_no_shared_peptides.fasta"

df_filtered = read_fasta_to_df(filename)


len(df_raw)
len(df_filtered)
len(df_raw)-len(df_filtered) # the difference between raw and without shared peptides.
len(df_filtered)/len(df_raw) # We filter away 26% of all proteins.

