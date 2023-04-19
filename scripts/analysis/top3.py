#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:19:23 2022

@author: ptruong
"""


import os 
import pandas as pd
import numpy as np
import scipy.stats as stats
from q_value import qvalues
import argparse


def top3(df):
    A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "1"].isna().sum(axis=1)<2]
    A = A.iloc[:,A.columns.get_level_values("sample_id") == "1"]
    B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "2"].isna().sum(axis=1)<2]
    B = B.iloc[:,B.columns.get_level_values("sample_id") == "2"]
    
    print("Finding all overlapping proteins in Sample A and Sample B.")
    # Find overlapping proteins
    overlapping_proteins = list(set(A.index) & set(B.index))
    A = A[A.index.isin(overlapping_proteins)]
    B = B[B.index.isin(overlapping_proteins)]
    
    print("Computing p-values using scipy.stats.ttest_ind() between sample A and B.")
    p_vals = stats.ttest_ind(A, B, axis = 1)[1]
    p_vals = pd.DataFrame(p_vals, columns = ["p"])
    print("Computing q-values.")
    p_vals["q"] = qvalues(p_vals)
    p_vals = pd.DataFrame(p_vals.values, index = A.index, columns = ["p", "q"])
    p_vals.sort_values("q",inplace =True)
    p_vals = p_vals.astype(float)
    
    print("Applying log2-transformation to sample sums")
    A = np.log2(A.sum(axis=1))
    B = np.log2(B.sum(axis=1))
    
    A.name = "1"
    B.name = "2"
    
    df_final = pd.concat([A, B, p_vals], axis = 1)
    print("Computing fold-change")
    df_final["log2(A,B)"] = df_final["1"] - df_final["2"]
    return df_final

parser = argparse.ArgumentParser(
    description='This script takes input from top3_formatter/diann.py and top3_formatter/osw.py and computes the protein quantity, p-values and q-values.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input', type=str,
                    help='Top3 input file. This file be generated using top3_formatter/diann.py and top3_formatter/osw.pyl')

parser.add_argument('--output', type=str, default = "top3_output.csv",
                    help='Name of the output file.')


# parse arguments from command line
args = parser.parse_args()
input_file = args.input
output = args.output

if __name__ == "__main__":
    print("Reading file: " + input_file)
    input_df = pd.read_csv(input_file, sep = "\t", index_col = [0,1], header = [0,1])
    print("Starting top3 computation...")
    res_df = top3(input_df)
    print("Printing output: " + output)
    res_df.to_csv(output, sep = "\t")
    print("Done!")
    















