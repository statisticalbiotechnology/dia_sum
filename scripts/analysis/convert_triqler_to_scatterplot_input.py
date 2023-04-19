#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:04:22 2022

@author: ptruong
"""


import pandas as pd
import numpy as np
from parsers.parse_triqler import parse_triqler
import argparse

def convert_triqler_to_scatterplot_input(triqler_results, triqler_input):
    """
    We need both triqler results file and triqler input data, because
    triqler results file does not contain absolute protein quantities
    (It only contains relative values).
    """
    df_triqler_res = parse_triqler(triqler_results)
    df_triqler_res = (df_triqler_res[~((df_triqler_res.iloc[:,df_triqler_res.columns.str.contains("Pedro")] == 1).sum(axis=1) == 6)]) # Remove proteins that are not updates at any sample

    df_triqler_input = pd.read_csv(triqler_input, sep = "\t")
    
    df_triqler_res["decoy"] = df_triqler_res.protein.map(lambda x:x.split("_")[0])
    df_triqler_res = df_triqler_res[df_triqler_res["decoy"] != "DECOY"] #dilter away decoy
    df_triqler_res["specie"] = df_triqler_res.protein.map(lambda x:x.split("_")[1])
    
    protein_list = df_triqler_res.protein.unique()
    
    df_triqler_input = df_triqler_input[df_triqler_input.proteins.isin(protein_list)]
    df_triqler_input["q"] = (1/np.e**df_triqler_input.searchScore) # reverse the searchScore to q-value
    
    means = []
    for i in df_triqler_input.condition.unique():
        means.append(df_triqler_input[df_triqler_input.condition == i].groupby("proteins").intensity.mean().rename(str(i)))
    
    df = pd.concat(means, axis = 1)
    df = np.log2(df)
    df["log2(A,B)"] = -df_triqler_res.set_index("protein").log2_fold_change
    df["specie"] = df.index.map(lambda x:x.split("_")[1])
    df["q"] = df_triqler_res.set_index("protein").q_value
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate triqler protein quantities and differential abundance results to scatter plot df.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--triqler_results', type=str,
                        help='input triqler result file.')
    parser.add_argument('--triqler_input', type = str,
                        help="input triqler input file.")
    parser.add_argument('--output', type=str, default = "triqler_scatter_format_input.csv",
                        help='output name.')

    # parse arguments from command line
    args = parser.parse_args()
    triqler_results = args.triqler_results
    triqler_input = args.triqler_input
    output = args.output
    
    df = convert_triqler_to_scatterplot_input(triqler_results, triqler_input)
    df.reset_index(inplace = True)
    df.rename({"proteins":"ProteinName"}, axis = 1, inplace = True)
    df = df[['specie', 'ProteinName', '1', '2', 'q', 'log2(A,B)']]
    df.to_csv(output, sep = "\t", index = False)
    

