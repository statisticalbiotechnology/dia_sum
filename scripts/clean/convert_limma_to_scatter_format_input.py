#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:31:45 2022

@author: ptruong
"""

import pandas as pd
import numpy as np 
import argparse

def convert_limma_to_scatter(limma_file):
    df_method = pd.read_csv(limma_file, sep = "\t")
    df_method = df_method.rename({"logFC":"log2(A,B)", "adj.P.Val":"q", "Unnamed: 0":"ProteinName"}, axis = 1)
    df_method["log2(A,B)"] = -df_method["log2(A,B)"]
    return df_method


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Concatenate limma results to scatter plot df.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--limma_results', type=str,
                        help='input limma result file.')
    parser.add_argument('--output', type=str, default = "limma_scatter_format_input.csv",
                        help='output name.')

    # parse arguments from command line
    args = parser.parse_args()
    limma_results = args.limma_results
    output = args.output
    
    df = convert_limma_to_scatter(limma_results)
    df.to_csv(output, sep = "\t", index = False)


