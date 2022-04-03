#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 21:50:17 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import argparse



def convert_diann_to_msqrob2(input_file, output, fdr_threshold = 0.01, qCol = "q"):
    print("Converting diann to msqrob2 format")
    print(f"fdr_threshold: {fdr_threshold}")
    df = pd.read_csv(input_file, sep = "\t")
    fdr_threshold = fdr_threshold
    df = df[df[qCol] < fdr_threshold]
    df["decoy"] = df["Protein.Ids"].str.contains("DECOY").copy(deep=True)
    df = df[df.decoy != True] # Filter away decoy as per msqrob2 manual
    df = df[~df.decoy.isna()] # Filter away decoy as per msqrob2 manual
    df = df.reset_index().drop("index", axis = 1)
    df = df.sort_values(by=qCol)
    df = df.groupby(["Stripped.Sequence", "Run"]).head(1) #127891
    
    df_pivot = df.pivot(index=["Stripped.Sequence","Protein.Ids"], columns = "Run", values = "Precursor.Quantity")
    df_pivot.reset_index(inplace=True)
    df_pivot["Proteins"] = df_pivot["Protein.Ids"]
    df_pivot.drop("Protein.Ids", axis = 1)
    df_pivot.to_csv(output, sep = "\t", index = False)
    print("convert_diann_to_msqrob2 Done!")
    

parser = argparse.ArgumentParser(
    description='Converts diann to MsqRob2 input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--fdr_threshold', type=float, default = 0.01,
                    help='fdr threshold to apply on qCol.')

parser.add_argument('--qCol', type=str, default = "q",
                    help='q value column to use.')

parser.add_argument('--output', type=str, default = "msqrob2_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
fdr_threshold = args.fdr_threshold
qCol = args.qCol
output = args.output

if __name__ == "__main__":
    convert_diann_to_msqrob2(input_file, output, fdr_threshold, qCol)


