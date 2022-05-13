#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:01:40 2022

@author: ptruong
"""


import pandas as pd
import numpy as np
import argparse

def convert_msstats_to_scatterplot(input_res, input_protein_quant):    
    df = pd.read_csv(input_protein_quant, sep = ",")
    df_ = pd.read_csv(input_res, sep = ",")
    df = df[["Protein", "LogIntensities", "GROUP"]]
    df = pd.concat([df_.set_index("Protein"), 
               pd.DataFrame(df[df.GROUP == 1].groupby("Protein").mean().LogIntensities).rename({"LogIntensities": "1"}, axis = 1),
               pd.DataFrame(df[df.GROUP == 2].groupby("Protein").mean().LogIntensities).rename({"LogIntensities": "2"}, axis = 1)],
              axis = 1)
    df.reset_index(inplace=True)
    df["specie"] = df.Protein.map(lambda x:x.split("_")[1])
    df = df[["specie", "Protein", "1", "2", "pvalue", "adj.pvalue", "log2FC"]].rename(
        {"Protein": "ProteinName", "log2FC":"log2(A,B)", "pvalue":"p", "adj.pvalue":"q"}, axis = 1)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate msstats protein quantities and differential abundance results to scatter plot df.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_res', type=str,
                        help='input result file.')

    parser.add_argument('--input_protein_quant', type=str,
                        help='input protein quant result file')

    parser.add_argument('--output', type=str, default = "msstats_scatter_format_input.csv",
                        help='output name.')

    # parse arguments from command line
    args = parser.parse_args()
    input_res = args.input_res
    input_protein_quant = args.input_protein_quant
    output = args.output
    
    df = convert_msstats_to_scatterplot(input_res, input_protein_quant)
    df.to_csv(output, sep = "\t", index = False)
    

