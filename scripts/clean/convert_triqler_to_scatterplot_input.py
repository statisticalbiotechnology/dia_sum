#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:03:20 2022

@author: ptruong
"""
import pandas as pd
import numpy as np
import argparse
from parsers.parse_triqler import parse_triqler

def convert_triqler_to_scatterplot(input_file):
    df = parse_triqler(input_file)
    df["2"] = df.iloc[:,df.columns.str.contains("2:")].mean(axis = 1)
    df["1"] = df.iloc[:,df.columns.str.contains("1:")].mean(axis = 1)
    df["decoy"] = df.protein.map(lambda x:x.split("_")[0])
    df = df[df["decoy"] != "DECOY"] #dilter away decoy
    df["specie"] = df.protein.map(lambda x:x.split("_")[1])
    df = df.rename({"protein":"ProteinName", "q_value":"q", "log2_fold_change":"log2(A,B)"}, axis = 1)
    df = df[["specie", "ProteinName", "1", "2", "q", "log2(A,B)"]]
    df.specie.unique()
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate triqler protein quantities and differential abundance results to scatter plot df.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_res', type=str,
                        help='input result file.')

    parser.add_argument('--output', type=str, default = "triqler_scatter_format_input.csv",
                        help='output name.')

    # parse arguments from command line
    args = parser.parse_args()
    input_res = args.input_res
    output = args.output
    df = convert_triqler_to_scatterplot(input_res)
    df.to_csv(output, sep = "\t", index = False)
    