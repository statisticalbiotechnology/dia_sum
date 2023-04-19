#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:12:41 2021

@author: ptruong
"""

import pandas as pd
import os 

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")


df = pd.read_csv("aligned.csv", sep = "\t")

split_agg_Frag = lambda x: len(x.split(";"))


df["aggr_Fragment_annotation_count"] = df.aggr_Fragment_Annotation.map(split_agg_Frag)

df["aggr_Fragment_annotation_count"].min()
df["aggr_Fragment_annotation_count"].max()
