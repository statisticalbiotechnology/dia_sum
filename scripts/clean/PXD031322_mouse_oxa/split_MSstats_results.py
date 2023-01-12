#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 15:49:36 2023

@author: ptruong
"""

import pandas as pd

df = pd.read_csv("LT_ST_Ctrl_msstats_results.csv", sep = ",")
df[df["Label"] == "LT-Ctrl"].to_csv("LT_Ctrl_msstats_results.csv", sep = ",", index = False)
df[df["Label"] == "ST-Ctrl"].to_csv("ST_Ctrl_msstats_results.csv", sep = ",", index = False)
df[df["Label"] == "ST-LT"].to_csv("ST_LT_msstats_results.csv", sep = ",", index = False)






