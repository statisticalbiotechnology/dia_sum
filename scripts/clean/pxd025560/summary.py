#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:24:47 2023

@author: ptruong
"""

import pandas as pd
import numpy as np
import re

df = pd.read_csv("pxd025560.report.tsv", sep = "\t")


df.Run.unique()
len(df.Run.unique())
sum(pd.DataFrame(df.Run.unique())[0].str.contains("5ug"))
sum(pd.DataFrame(df.Run.unique())[0].str.contains("reprep"))
sum(pd.DataFrame(df.Run.unique())[0].str.contains("DIA_HF"))

# Early-stage cohort - tumour samples from 192 patients with lung cancer.
# Late-stage cohort - inoperable samples from 84 samples.


# Fig 8
meta = pd.read_csv("metadata.csv", sep = "/")
meta.Histology.unique()
mapper = meta.set_index("SampleID")["Histology"].to_dict()

df["sample_id"] =  df.Run.map(lambda x: re.findall('\d+',x.split("_")[-1])[0]).astype(int)
df["sample_class"] = df.sample_id.map(mapper)
df.sample_class



