#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 13:28:18 2022

@author: ptruong
"""

import os
import pandas as pd
import numpy as np

for i in os.listdir():
    print(i)
    df = pd.read_csv(i, sep = "\t", usecols = ["decoy", "ProteinName", "FullPeptideName", "filename", "m_score"])
    df = df[df.decoy != 1]
    #df = df[df.m_score < 0.000794328234724281] # m_score cutoff with target FDR 0.01
    n_proteins = len(df.ProteinName.unique())
    n_peptides = len(df.FullPeptideName.unique())
    filename = df.filename.unique()[0]
    print(filename)
    print("proteins:" + str(n_proteins))
    print("peptides:" + str(n_peptides))
    print("----")






