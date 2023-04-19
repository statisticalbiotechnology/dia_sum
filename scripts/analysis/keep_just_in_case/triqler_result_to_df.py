#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 19:30:05 2021

@author: ptruong
"""

import pandas as pd
import numpy as np

f = open("proteins.tsv")
cols = f.readline().split("\n")[0].split("\t")

val_array = []
for line in f:
    line = line.split("\n")[0].split("\t")
    vals = line[:len(cols)-1]
    peptides = line[len(cols)-1:]
    peptides = ";".join(peptides)
    vals.append(peptides)
    val_array.append(vals)

df = pd.DataFrame(val_array, columns = cols)

df.to_csv("proteins_proc.csv", sep = "\t")
df= pd.read_csv("proteins_proc.csv", sep = "\t")


(df["q_value"]  < 0.05).sum() # 966 differentially abundant protein






