#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 12:16:35 2021

@author: ptruong
"""

import os 
import re
import pandas as pd
import numpy as np

# read in .fasta file and count shared peptides

os.chdir("/home/ptruong/git/dia_sum")
filename = "database/2021-06-07/UP00000625_UP000002311_UP000005640.fasta"
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


counted_df

counted_df[counted_df.protein > 1]


counted_df[counted_df.protein == 1]

# Conclusion

print(f"15167 sequences are shared in > 1 Proteins. ({round(len(counted_df[counted_df.protein > 1])/ len(counted_df), 6)} of proteins. ))")
print(f"682301 sequences has unique proteins.({round(len(counted_df[counted_df.protein == 1])/ len(counted_df), 6)} of proteins.))")
    
df
    
seqeunces_in_more_than_one_protein = counted_df[counted_df > 1].dropna().index
df_shared_peptides_protein = df[df.sequence.isin(seqeunces_in_more_than_one_protein)]


seqeunces_one_protein = counted_df[counted_df.protein == 1].dropna().index
df_non_shared_peptides_protein = df[df.sequence.isin(seqeunces_one_protein)]


seqeunces_in_more_than_one_protein.isin(seqeunces_one_protein).sum()


len(df_shared_peptides_protein)/len(df)

# Conclusion

print(f"{len(counted_df[counted_df.protein > 1])} sequences are shared in > 1 Proteins. ({round(len(counted_df[counted_df.protein > 1])/ len(counted_df), 6)} of proteins. ))")
print(f"{len(counted_df[counted_df.protein == 1])} sequences has unique proteins.({round(len(counted_df[counted_df.protein == 1])/ len(counted_df), 6)} of proteins.))")

# Conclusion
print(f"{len(df.protein.unique())} proteins in .FASTA")
print(f"{len(df_shared_peptides_protein.protein.unique())} proteins have shared peptides. ({round(len(df_shared_peptides_protein.protein.unique())/len(df.protein.unique()), 6)} of proteins. ))")
#print(f"{len(df_non_shared_peptides_protein)} proteins have unique peptides.({round(len(df_non_shared_peptides_protein.protein.unique()) / len(df.protein.unique()), 6)} of proteins.))")

    
df_non_shared_peptides_protein[df_non_shared_peptides_protein.protein.isin(df_shared_peptides_protein.protein)]


shared = df_shared_peptides_protein
non_shared = df_non_shared_peptides_protein´



######


shared_resampled = shared.groupby("protein").sample(1)
counts = shared_resampled.groupby("protein").count()
while counts.sequence.max() > 1:
    shared_resampled = shared.groupby("pro").sample(1)
    counts = shared_resampled.groupby("protein").count()


count = shared.groupby("sequence")["protein"].count()
count[count.index == "NIALIFEK"]

resample = shared.groupby("sequence").sample(1)
count = resample.groupby("sequence")["protein"].count()
shared_resampled = shared[shared.protein.isin(resample.protein)]
while shared_resampled.groupby("sequence").count().protein.max() > 1:
    print(shared_resampled.groupby("sequence").count().protein.max())
    resample = shared_resampled.groupby("sequence").sample(1)
    count = resample.groupby("sequence")["protein"].count()
    shared_resampled = shared_resampled[shared_resampled.protein.isin(resample.protein)]
# There seem to be some protein which always have overlapping peptides... perhaps discard these?

shared_resampled.to_csv("shared_resampled.csv", sep = "\t", index = False)


shared_resampled_ = pd.read_csv("shared_resampled.csv", sep = "\t")


shared_count = shared_resampled.groupby("sequence").count().protein

single_protein_sequence = shared_count[shared_count < 2].index # dropp sequences whose peptides match to more than 1 protein.

shared_resampled[shared_resampled.sequence.isin(single_protein_sequence)] #8867
shared_resampled #25994

test = shared.groupby("sequence").sample(1)
test_counts = shared_resampled.groupby("sequence").count().protein
test_counts[test_counts<2]
len(single_protein_sequence)


unique_protein_sequences = list(single_protein_sequence) + list(df_non_shared_peptides_protein.sequence)

unique_sequence_protein =