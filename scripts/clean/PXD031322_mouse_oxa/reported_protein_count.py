#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 14:41:52 2022

@author: ptruong
"""

import pandas as pd
import numpy as np 



df = pd.read_excel("reported_6_subclusters.xlsx", header = 1)
df = df[~df.Condition.isin(["Others"])]

pd.DataFrame(df.groupby("Condition").count()["Protein.Ids"].values, 
             index = df.groupby("Condition").count()["Protein.Ids"].index.values,
             columns = ["Reported"]).to_csv("reported_protein_count.tsv", sep = "\t")



