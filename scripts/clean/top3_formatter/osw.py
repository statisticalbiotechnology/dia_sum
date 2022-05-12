#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 16:53:24 2022

@author: ptruong
"""

import os
import pandas as pd
import numpy as np
import argparse

experiment_id_mapper = lambda x: x.split("_")[5]
sample_id_mapper = lambda x: x.split("_")[8] 
specie_mapper = lambda x: x.split("_")[-1]

# OSW
def avg_3_largest_precursor_on_run_level(filename, m_score_treshold):  
    print("Reading in: " + filename)
    df = pd.read_csv(filename, sep = "\t")
    df = df[df.decoy != 1]
    df = df[df.m_score < m_score_treshold] 
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    sample_id = df.sample_id.unique()[0]
    experiment_id = df.experiment_id.unique()[0]     
    def top3(df):
        df = (df.groupby('ProteinName')['Intensity'].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 2 else np.nan)
                  .reset_index())
        return df
    df_reduced = df[["ProteinName", "Intensity"]]
    df_protein = top3(df_reduced)
    df = df_protein
    df["specie"] = df.ProteinName.map(specie_mapper)
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df = df.set_index(["specie", "ProteinName"])
    df = pd.DataFrame(df.values, columns = midx, index = df.index)
    
    return df

def avg_3_largest_precursor(result_file_directory, m_score_threshold_file):
    m_score_threshold = float(open(m_score_threshold_file).read())

    dfs = []
    for file in os.listdir(result_file_directory):
        dfs.append(avg_3_largest_precursor_on_run_level(result_file_directory + "/" + file, m_score_treshold=m_score_threshold))
    df = pd.concat(dfs, axis = 1)
    return df

parser = argparse.ArgumentParser(
    description='This scripts formats OpenSwath output with dscore to take the avg of three largest PSM as peptide quantity. The Openswath outputs each sample as a file, therefore we input a file directory to this script. Note: The OpenSwath output with dscore need to be in a seperate directory for themselves.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--file_dir', type=str,
                    help='OpenSwath output with dscore file directory. Note: The OpenSwath output with dscore need to be in a seperate directory for themselves.')

parser.add_argument('--m_score_threshold_file', type=str,
                    help='m_score level to cut-off at. Note: use mscore4pepfdr() in Bioconductor SWATH2stats to compute m_score for desired peptide level fdr.')

parser.add_argument('--output', type=str, default = "osw_top3_input_formatted.csv",
                    help='Name of the output file.')


# parse arguments from command line
args = parser.parse_args()
file_dir = args.file_dir
m_score_threshold_file = args.m_score_threshold_file
output = args.output

if __name__ == "__main__":
    df = avg_3_largest_precursor(file_dir, m_score_threshold_file = m_score_threshold_file)
    print("Generating output file: " + output)
    df.to_csv(output, sep = "\t")