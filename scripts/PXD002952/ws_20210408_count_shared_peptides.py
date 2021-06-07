#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 18:59:43 2021

@author: ptruong
"""

import re

# read in .fasta file and count shared peptides

filename = "UP00000625_UP000002311_UP000005640.fasta"

file = open(filename, "r")


#line = file.readline()
protein_list = []
sequence_list = []
#protein_dict = {}
for line in file: 
    if line[0] == ">":
        protein = line # line.split("_")[0].split("|")[1] + "_" + line.split("_")[1].split(" ")[0]
        #protein_list.append(protein)
        pass
    else:
        sequence = line.rstrip()
        split_sequence = re.split(r"(?<=[KR])", sequence)
        split_sequence = list(dict.fromkeys(split_sequence))
        #protein_dict[protein] = split_sequence
        sequence_list += split_sequence
        protein_list += [protein for i in range(len(split_sequence))]
        

import time
import pandas as pd
import numpy as np
start = time.time()
df = pd.DataFrame(np.array([protein_list, sequence_list]).T, columns = ["protein", "sequence"])
end = time.time()
print(end -start)
df["seq_length"] = df.sequence.str.len()
df = df[df["seq_length"] > 7]
df.drop("seq_length", axis = 1, inplace = True)

counted_df = df.groupby("sequence").count().sort_values(by = "protein", ascending = False)

counted_df
counted_df[counted_df.protein > 1]

#df.drop_duplicates()
#df_dropped = df.drop_duplicates()

#counted_df_dropped = df_dropped.groupby("sequence").count().sort_values(by = "protein", ascending = False)

#counted_df_dropped
#counted_df[counted_df.protein > 1]
