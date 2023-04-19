#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 22:21:49 2022

@author: ptruong
"""
import os 
import numpy as np
import pandas as pd
import argparse


def convert_osw_to_triqler(input_file, output):
    df = pd.read_csv(input_file, sep = "\t")
    # filename has different formatting, we need to change number or implement regex.
    experiment_id_mapper = lambda x: x.split("_")[5]
    sample_id_mapper = lambda x: x.split("_")[8] #hye124 
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    df_triq = df[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]
    df_triq = df_triq.rename(columns={"experiment_id": "run", "sample_id": "condition", "Charge": "charge", 
                            "m_score":"searchScore", "Intensity":"intensity", "FullPeptideName":"peptide",
                            "ProteinName":"proteins"}, errors="raise")
    df_triq.searchScore += 0.000001
    df_triq["searchScore"] = -np.log10(df_triq["searchScore"])
    df_triq.to_csv(output, sep = "\t", index=False)

        
 
parser = argparse.ArgumentParser(
    description='Convert OSW to Triqler input format.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--input_file', type=str,
                    help='input file name.')

parser.add_argument('--output', type=str, default = "osw_triqler_input.csv",
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
input_file = args.input_file
output = args.output
      
if __name__ == "__main__":
    convert_osw_to_triqler(input_file, output)



