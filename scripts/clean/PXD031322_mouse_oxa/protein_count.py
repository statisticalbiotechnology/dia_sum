#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 14:47:58 2022

@author: ptruong
"""

import pandas as pd
import numpy as np

reported_file = "reported_protein_count.tsv"
triqler_file = "triqler_protein_count_fc_0.415.tsv"
top3_file = "top3_protein_count.tsv"
msstats_file = "msstats_protein_count.tsv"
msqrob2_file = "msqrob2_protein_count.tsv"

reported = pd.read_csv(reported_file, sep = "\t").rename({"Unnamed: 0": "Condition"}, axis = 1).set_index("Condition")
triqler = pd.read_csv(triqler_file, sep = "\t").rename({"Unnamed: 0": "Condition"}, axis = 1).set_index("Condition")
top3 = pd.read_csv(top3_file, sep = "\t").rename({"Unnamed: 0": "Condition"}, axis = 1).set_index("Condition")
msstats = pd.read_csv(msstats_file, sep = "\t").rename({"Unnamed: 0": "Condition"}, axis = 1).set_index("Condition")
msqrob2 = pd.read_csv(msqrob2_file, sep = "\t").rename({"Unnamed: 0": "Condition"}, axis = 1).set_index("Condition")

protein_count_table = pd.concat([reported, triqler, top3, msstats, msqrob2], axis = 1)

protein_count_table.to_csv("protein_count_triqler_fc_0.415.tsv", sep = "\t")


