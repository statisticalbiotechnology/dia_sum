#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 15:07:47 2021

@author: ptruong
"""



import pandas as pd
import numpy as np
import os 

m_score_treshold = 1.00
directory = "/hdd_14T/data/PXD002952/20210805_osw_run"
output_name = "triqler_input_filtered_mscore_0.01.csv"

def osw_output_to_triqler(directory, output_name, m_score_treshold = 1.00):
    df_all = pd.DataFrame()
    for file in os.listdir(directory):
        if file[-10:] == "dscore.csv":
            df = pd.read_csv(directory + "/" + file, sep = "\t")
            print(file)
            print(len(df))
            df_all = pd.concat([df_all, df], axis = 0)
    df_all = df_all.reset_index().drop("index", axis = 1)
    df_all = df_all[df_all.m_score < m_score_treshold] #tresholding for msqrobsum run.
    
    # filename has different formatting, we need to change number or implement regex.
    experiment_id_mapper = lambda x: x.split("_")[5]
    sample_id_mapper = lambda x: x.split("_")[8] #hye124 
    df_all["experiment_id"] = df_all["filename"].map(experiment_id_mapper)
    df_all["sample_id"] = df_all["filename"].map(sample_id_mapper)
    
    df_triq = df_all[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]
    df_triq = df_triq.rename(columns={"experiment_id": "run", "sample_id": "condition", "Charge": "charge", 
                            "m_score":"searchScore", "Intensity":"intensity", "FullPeptideName":"peptide",
                            "ProteinName":"proteins"}, errors="raise")
    df_triq["searchScore"] = -np.log10(df_triq["searchScore"])
    df_triq.to_csv(directory + "/" + output_name, sep = "\t", index=False)
    
if __name__ == "__main__":
    osw_output_to_triqler(directory, output_name, m_score_treshold = 0.01)