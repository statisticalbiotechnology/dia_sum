#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 20:17:16 2021

@author: ptruong
"""

import pandas as pd 
import numpy as np

aligned = pd.read_csv("aligned.csv", sep = "\t")

species_map = lambda x : x.split("_")[-1]
run_map = lambda x : x.split("_")[5]
sample_map = lambda x : x.split("_")[8]


species_list = aligned.ProteinName.map(species_map).unique()
run_list = aligned.filename.map(run_map).unique()
sample_list = aligned.filename.map(sample_map).unique()

aligned["specie"] = aligned.ProteinName.map(species_map)
aligned["run"] = aligned.filename.map(run_map)
aligned["sample"] = aligned.filename.map(sample_map)

# treshold

midx = pd.MultiIndex.from_tuples(tuples, names=["specie", "run", "sample"])


tuples = list(aligned[["specie", "run", "sample"]].T.values)
