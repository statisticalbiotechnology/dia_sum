#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 10:22:38 2022

@author: ptruong

http://localhost:8888/notebooks/result/2021-09-13_library_peptide_filtering/20210913_filter_fasta_and_analysis_of_databse.ipynb

"""

import os 
import re
import pandas as pd
import numpy as np


os.chdir("/home/ptruong/git/dia_sum")

#human_fasta = "2021-06-15-decoys-reviewed-contam-UP000000625.fas"
#yeast_fasta = "2021-06-15-decoys-reviewed-contam-UP000002311"
#ecoli_fasta = "2021-06-15-decoys-reviewed-contam-UP000005640"

#filename = "database/2021-06-15/" + human_fasta
#filename = "database/napedro_3mixed_human_yeast_ecoli_20140403_iRT_reverse.fasta"
#filename = "database/2021-06-07/UP00000625_UP000002311_UP000005640.fasta"
#filename = "database/2021-09-15_malaria_yeast_PXD017705/UP000002311_UP000001450.fasta"
filename = "database/2022-05-25_pxd028185/SupplementaryFiles files/libraries and fastas/uniprot_human_25apr2019_yeastENO1.fasta"
# filename = "/home/ptruong/git/dia_sum/database/2021-11-02_PXD026600_ecoli/UP000000625_reviewed_ups1-ups2-sequences_2021-11-02.fasta"


def get_no_isoforms_proteins(filename):
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
    
    
    
    mofidy_IL = True #Set the variable
    if mofidy_IL == True:
        modify_IL_function = lambda x: x.replace("I", "L")
        df["sequence"] = df.sequence.map(modify_IL_function)
        print("All I and L are now L")
    else:
        print("No I/L modification")
        
    def decoy_map(protein):
        if protein.split("_")[0] == ">reverse":
            return True
        else:
            return False
        
        
    df["decoy"] = df.protein.map(decoy_map)
    
    df = df[df.decoy == False]
    specie_map = lambda x: x.split("_")[1].split(" ")[0]
    df["specie"] = df.protein.map(specie_map)
    
    n_protein = len(df.protein.unique())
    n_protein_ecoli  = len(df[df["specie"] == "ECOLI"].protein.unique())
    n_protein_human = len(df[df["specie"] == "HUMAN"].protein.unique())
    n_protein_yeast = len(df[df["specie"] == "YEAST"].protein.unique())
    
    print("Unfiltered protein statistics:")
    print(f"All proteins: {n_protein}")
    print(f"ECOLI proteins: {n_protein_ecoli}")
    print(f"HUMAN proteins: {n_protein_human}")
    print(f"YEAST protein: {n_protein_yeast}")
    
    df["seq_length"] = df.sequence.str.len()
    df = df[df["seq_length"] > 7]
    df.drop("seq_length", axis = 1, inplace = True)
    df.drop_duplicates(inplace=True)
    
    counted_df = df.groupby("sequence").count().sort_values(by = "protein", ascending = False)
    single_protein_seq = counted_df[counted_df.protein == 1]
    df_single_seq = df[df.sequence.isin(single_protein_seq.index)]
    keep_list = df_single_seq.protein.unique()
    return keep_list    
    

def remove_shared_peptide_proteins(fasta_input, fasta_output):
    
    keep_list = get_no_isoforms_proteins(filename = fasta_input)
    
    
    for line in fasta_input:
        if line[0] == ">":
            protein = line 
            fasta_output.write(protein)
        else:
            if protein in keep_list:
                sequence = line
                fasta_output.write(sequence)
