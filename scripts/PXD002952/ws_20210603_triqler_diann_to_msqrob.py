#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:14:25 2021

@author: ptruong
"""


import os 
# git 
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_to_msqrobsum_converter import get_mSqRobSum_input
# MSFragger 
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")

triqler = pd.read_csv("triqler_input_diann.csv", sep = "\t")
expr_, fd_, pd_ = get_mSqRobSum_input(triqler)
os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger/msqrob_input")
remap_sample = lambda s : "X"+".".join(s.split("-"))
pd_["sample"] = pd_["sample"].map(remap_sample)

expr_.to_csv("expr.csv", sep = "\t", index = False)
fd_.to_csv("fd.csv", sep = "\t", index = False)
pd_.to_csv("pd.csv", sep = "\t", index = False)
