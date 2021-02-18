#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:39:00 2021

@author: ptruong
"""

import os 
import pandas as pd

df = pd.read_csv("hroest_K120808_Strep0%PlasmaBiolRepl1_R01_SW.mzML.gz.tsv", sep = "\t")

# Merge all .tsv
df = pd.DataFrame()
for i in os.listdir():
    if i[-4:] == ".tsv":
        df = df.append(pd.read_csv(i, sep = "\t"))
        
df.to_csv("merged.tsv", sep = "\t", index=False)






