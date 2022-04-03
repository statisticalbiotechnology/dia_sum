#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 07:49:33 2022

@author: ptruong
"""


import numpy as np
import pandas as pd
import argparse

def diann_to_triqler(filename, qvalue_treshold = 1.00, qCol = "q"):
    df = pd.read_csv(filename, sep = "\t")
    
    run_mapper = lambda x : x.split("_")[5]
    condition_mapper = lambda x : x.split("_")[8]
    
    df["run"] = df["Run"].map(run_mapper)
    df["condition"] = df["Run"].map(condition_mapper)
    df["charge"] = df["Precursor.Charge"]
    #df["searchScore"] = df["CScore"]
    df["searchScore"] = df[qCol]
    df["intensity"] = df['Precursor.Quantity']
    df["peptide"] = df["Stripped.Sequence"]
    df["proteins"] = df["Protein.Ids"]
    
    df_triq = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    return df_triq

def main(input_file, output, qCol = "q"):
    df_triq = diann_to_triqler(input_file)
    df_triq.intensity = np.log(df_triq.intensity) #if required to log intensity.
    df_triq.to_csv(output, sep = "\t", index = False)

parser = argparse.ArgumentParser(
    description='Converts diann to Triqler input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--qCol', type=str, default = "q",
                    help='q value column to use.')

parser.add_argument('--output', type=str, default = "msqrob2_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
qCol = args.qCol
output = args.output

if __name__ == "__main__":
    main(input_file, output, qCol = "q")

