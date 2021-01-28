#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 12:00:15 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

def get_log2FC_matrix(df, sample1, sample2, specie):
    """
    get log2FC matrix for all runs ans proteins for a method (triq, spec) and between
    samples, for a specie.
    
    df - melted triq or spec from extract_from_melt.py
    sample - sample
    specie - specie
    """
    proteins = get_proteins(df, specie) 
    runs = ["R0"+str(i) for i in range(1,6)]
    log2FC_array = []
    start=time.time()
    i = 0
    for protein in proteins:
        if i%10 == 0:
            print(str(i) + "/" + str(len(proteins)))
        log2FC_array.append(parallel_get_protein(df, protein, sample1, sample2, specie))
        i+=1
    end=time.time()
    print(end-start)
    df_log2FC = pd.DataFrame(log2FC_array, index=proteins, columns=runs)
    return df_log2FC



