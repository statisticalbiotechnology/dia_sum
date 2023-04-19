#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 01:04:22 2022

@author: ptruong
"""

import pandas as pd
import numpy as np


df = pd.read_csv("lib_converted_decoy.tsv", sep = "\t")
df["decoy_map_col"] = df["ProteinId"].map(lambda x:x.split("_")[0]) == "DECOY"

def map_decoy(x):
    if x == True:
        return "DECOY_"
    else:
        return ""

df["decoy_map_col"] = df["decoy_map_col"].map(map_decoy)
df["UniprotId"] = df["decoy_map_col"] + df["UniprotId"]
df.drop("decoy_map_col", axis = 1, inplace = True)

df.to_csv("lib_converted_decoy_added_tag.tsv", sep = "\t", index = False)








