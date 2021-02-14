#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:09:17 2021

@author: ptruong
"""


import pandas as pd 


df = pd.read_csv("aligned_osw.tsv", sep = "\t")

triq_vals = df[["filename", "run_id", "Charge", "m_score", "Intensity", "Sequence", "ProteinName"]].values

col_names = ["run","condition","charge",	"searchScore", "intensity",	"peptide",	"proteins"]

triq = pd.DataFrame(triq_vals, columns = col_names)

triq.to_csv("triqler_formatted.csv", sep = "\t", index = False)
