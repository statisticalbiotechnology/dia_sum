#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 13:38:46 2022

@author: ptruong
"""


import numpy as np
import pandas as pd

def fdr(p_vals):
    "BH-corretion"
    "https://stackoverflow.com/questions/25185205/calculating-adjusted-p-values-in-python"
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr




df = pd.read_csv("report.pg_matrix.tsv", sep = "\t")
p = pd.read_csv("limma_p_values.csv", sep = "\t")
df["p"] = p
df["q"] = fdr(p.values.flatten())
df.to_csv("report.pg_matrix.p.tsv", sep = "\t")

