#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 17:35:50 2022

@author: ptruong
"""

import pandas as pd
import numpy as np

reported =  pd.read_excel("1-s2.0-S1874391922002068-mmc2.xlsx", header = 1)
reported = reported[reported.Condition != "Not_sig"]

triqler_file = "triqler_input_report_0.5.tsv"


def parse_triqler(triqler_output_file):
    """
    Parses triqler output format to pandas dataframe.
    """
    f = open(triqler_output_file, "r")
    lines = f.readlines()
    line = lines.pop(0)
    cols = line.split("\n")[0].split("\t")[:]
    n_cols = len(cols)
    
    data_array = []
    for line in lines:
        line = line.split("\n")[0].split("\t")
        vals = line[:n_cols-1]
        peptides = ";".join(line[n_cols-1:])
        data = vals + [peptides]
        data_array.append(data)
    df = pd.DataFrame(data_array, columns = cols)

    df = pd.concat([df[["protein", "peptides"]], df.drop(["protein", "peptides"], axis = 1).astype(float)], axis = 1)
    
    return df


triqler_file = "proteins.tsv"

triq = parse_triqler(triqler_file)
triq = triq[triq.q_value < 0.05]

triq["protein_name"] = triq.protein.map(lambda x:x.split("_")[0])


intersection = list(set(triq.protein_name) & set(reported.Gene_Symbol.str.upper()))

len(intersection)

len(triq.protein_name)
len(reported.Gene_Symbol.str.upper())
triq.protein_name.to_csv("triq_protein_name.tsv")
reported.Gene_Symbol.str.upper().to_csv("reported_protein_name.tsv")
# Triqler has gene name
# reported has uniprot prot name


