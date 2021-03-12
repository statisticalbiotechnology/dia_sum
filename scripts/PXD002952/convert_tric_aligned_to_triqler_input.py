#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 22:24:57 2021

@author: ptruong
"""


import pandas as pd
import numpy as np

df = pd.read_csv("aligned.csv", sep = "\t")

experiment_id_mapper = lambda x: x.split("_")[5]
#sample_id_mapper = lambda x: x.split("_")[7]
sample_id_mapper = lambda x: x.split("_")[9] #hye124 
df["experiment_id"] = df["filename"].map(experiment_id_mapper)
df["sample_id"] = df["filename"].map(sample_id_mapper)
#df["Intensity"] = np.log2(df["Intensity"])

df_triq = df[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]

df_triq = df_triq.rename(columns={"experiment_id": "run", "sample_id": "condition", "Charge": "charge", 
                        "m_score":"searchScore", "Intensity":"intensity", "FullPeptideName":"peptide",
                        "ProteinName":"proteins"}, errors="raise")
df_triq.to_csv("triqler_input.csv", sep = "\t", index=False)








