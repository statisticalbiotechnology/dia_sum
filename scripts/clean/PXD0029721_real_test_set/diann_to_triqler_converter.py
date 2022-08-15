#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 19:50:40 2022

@author: ptruong
"""

import numpy as np
import pandas as pd

def convert_diann_to_triqler(filename, output):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "Q.Value",
                                "Precursor.Quantity", "Stripped.Sequence",
                                "Protein.Ids"])
    df["condition"] = df.Run.map(lambda x: x.split("_")[-2])
    df.rename({"Run":"run", "Precursor.Charge":"charge", "Q.Value":"searchScore",
               "Precursor.Quantity":"intensity", "Stripped.Sequence":"peptide",
               "Protein.Ids":"proteins"}, axis = 1, inplace = True)
    df = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    df["searchScore"] = -np.log(df["searchScore"])
    df = df.dropna()
    df.to_csv(output, sep = "\t", index = False)


convert_diann_to_triqler(filename = "report.tsv", 
                         output = "triqler_input_report.tsv")







