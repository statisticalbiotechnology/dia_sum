#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 08:06:14 2022

@author: ptruong
"""


import os 
os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")

import pandas as pd
from time import time

df_example = pd.read_csv("example_disaggregate_output_format.csv", sep = ",")

# count sequences
df_example.groupby("PeptideSequence").count()

# count fragment ions
df_example.groupby("FragmentIon").count()

# Investigate a specific fragment ion
df_example[df_example["FragmentIon"] == 1154]


# A lot of the data consist of very low or 0 intensity fragment ions.
# These data points can most likely be filtered away.

