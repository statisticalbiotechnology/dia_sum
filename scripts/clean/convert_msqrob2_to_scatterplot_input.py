#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:02:26 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import argparse


def convert_msqrob2_to_scatterplot(input_res, input_protein_quant):
    df = pd.read_csv(input_protein_quant, sep = ",").rename({'Unnamed: 0':"Proteins"}, axis = 1).set_index("Proteins")
    df_ = pd.read_csv(input_res, sep = ",").rename({'Unnamed: 0':"Proteins"}, axis = 1).set_index("Proteins")
    df = pd.concat([df ,df_], axis=1, join='inner')
    df["2"] = df.iloc[:,df.columns.str.contains("Sample_2")].mean(axis = 1)
    df["1"] = df.iloc[:,df.columns.str.contains("Sample_1")].mean(axis = 1)
    
    df = df[["1", "2", "pval", "adjPval", "logFC"]].reset_index()
    df["specie"] = df.Proteins.map(lambda x:x.split("_")[1])
    df = df[["specie", "Proteins", "1", "2", "pval", "adjPval", "logFC"]].rename({"Proteins": "ProteinName", 
                                                                                  "logFC":"log2(A,B)",
                                                                                  "pval":"p",
                                                                                  "adjPval":"q"}, axis = 1)
    return df


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='Concatenate msqrob2 protein quantities and differential abundance results to scatter plot df.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_res', type=str,
                        help='input result file.')

    parser.add_argument('--input_protein_quant', type=str,
                        help='input protein quant result file')

    parser.add_argument('--output', type=str, default = "msqrob2_scatter_format_input.csv",
                        help='output name.')

    # parse arguments from command line
    args = parser.parse_args()
    input_res = args.input_res
    input_protein_quant = args.input_protein_quant
    output = args.output


    df = convert_msqrob2_to_scatterplot(input_res, input_protein_quant)
    df.to_csv(output, sep = "\t", index = False)
    
    
    
    
