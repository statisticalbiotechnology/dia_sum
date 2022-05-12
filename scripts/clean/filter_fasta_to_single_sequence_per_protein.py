#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:10:02 2022

@author: ptruong
"""

import re
import pandas as pd
import numpy as np
import argparse


def filter_fasta_to_single_sequence_per_protein(filename, output):
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
    single_protein_seq = counted_df[counted_df.protein == 1]
    df_single_seq = df[df.sequence.isin(single_protein_seq.index)]
    keep_list = df_single_seq.protein.unique()
    
    file = open(filename, "r")
    
    file_w = open(output, "w")
    
    for line in file:
        if line[0] == ">":
            protein = line 
            file_w.write(protein)
        else:
            if protein in keep_list:
                sequence = line
                file_w.write(sequence)
            
parser = argparse.ArgumentParser(
    description='Formats the .fasta file to single sequence per protein fasta.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--filename', type=str,
                    help='input fasta name.')

parser.add_argument('--output', type=str,
                    help='output fasta name.')

# parse arguments from command line
args = parser.parse_args()
filename = args.filename
output = args.output

if __name__ == "__main__":
    filter_fasta_to_single_sequence_per_protein(filename, output)

   
    
    
    
