#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 19:21:14 2021

@author: ptruong
"""

import os 
import pandas as pd
import numpy as np
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

from triqler_to_msqrobsum_converter import get_mSqRobSum_input


os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")



triqler = pd.read_csv("triqler_input.csv", sep = "\t")
expr_, fd_, pd_ = get_mSqRobSum_input(triqler)
# apply log2 to expr_ because mSqRobSum want log2 input.
#expr_ = pd.concat([expr_["peptide"], expr_.select_dtypes(include=['float64']).apply(np.log2)], axis = 1)
sample_to_r_mapper = lambda x : "X" + ".".join(x.split("-"))
pd_["sample"] = pd_["sample"].map(sample_to_r_mapper)
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/mSqRobSum_input")

expr_.to_csv("expr.csv", sep = "\t", index = False)
fd_.to_csv("fd.csv", sep = "\t", index=False)
pd_.to_csv("pd.csv", sep = "\t", index=False)
#pd.DataFrame(expr_.columns[1:]).T.to_csv("expr_col.csv", index = False, header = False, sep = "\t")



