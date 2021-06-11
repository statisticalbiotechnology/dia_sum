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


df_shared_peptides_protein


df_shared_peptides_protein.protein.unique()[0]

prot = shared[shared.sequence == "NIALIFEK"].iloc[0,:].protein

df[df.protein == prot]

df[df.protein == prot].sequence

