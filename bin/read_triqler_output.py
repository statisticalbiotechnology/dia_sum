#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 14:43:35 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

def read_triqler_protein_output_to_df(filename = "proteins.1vs2.tsv"):
    """
    Function is needed because protein.<x>vs<y>.tsv has tab-seperated peptides, 
    makes cases the read-in module in pandas to bug.
    """

    f = open(filename, "r")
    
    header = f.readline().split("\n")[0].split("\t")
    
    len_header = len(header)
    
    lines = []
    for line in f:
        line = line.split("\n")[0].split("\t")
        vals = line[:len_header - 1]
        peptides = line[(len_header - 1 ):]
        peptides = ";".join(peptides)
        vals.append(peptides)
        lines.append(vals)
    
    df = pd.DataFrame(lines, columns = header)    
    df[header[:2] + header[3:-1]] = df[header[:2] + header[3:-1]].apply(pd.to_numeric)

    return df




