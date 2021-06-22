#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 16:43:42 2021

@author: ptruong
"""

import pandas as pd
import numpy as np



def df_batch(filename, n_batches):
    df = pd.read_csv(filename)

    batches = np.array_split(df, n_batches)
    
    for batch in range(n_batches):
        batches[batch].to_csv(filename[:-4]+ f"_batch_{batch+1}" + ".csv", sep = ",", index = False)

filename = "no_shared_UP000005640_human.fasta.trypsin.z3_nce33.csv"
df_batch(filename, n_batches = 4)
    