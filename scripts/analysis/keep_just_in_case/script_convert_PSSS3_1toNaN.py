#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 13:06:04 2020

@author: ptruong
"""

import pandas as pd
import numpy as np 
df = pd.read_csv("PSSS3_triqlerFormatted_nonShared.csv", sep="\t")
df = df.replace(1, np.nan)
df.to_csv("PSSS3_triqlerFormatted_nonShared_1toNaN.csv", sep="\t", index=False)

