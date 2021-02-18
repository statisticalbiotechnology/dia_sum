#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:23:31 2021

@author: ptruong
"""

import pandas as pd

df = pd.read_csv("merged.tsv", sep = "\t")

len(df.run_id.unique())

x = df.filename.unique()[0]

def get_run(x):
    file_name_components = x.split("/")[-1].split("_")
    technical_replicate = file_name_components[-2]
    biological_replicate = file_name_components[-3].split("%")[1].split("PlasmaBiol")[-1]
    condition = file_name_components[-3].split("%")[0].split("Strep")[-1]

    run = "BR_" + biological_replicate + "_" + "TR" + "_" + technical_replicate + "_" + str(condition)
    return run

def get_condition(x):
    file_name_components = x.split("/")[-1].split("_")
    condition = file_name_components[-3].split("%")[0].split("Strep")[-1]
    if condition == "0":
        condition = "A"
    elif condition == "10":
        condition = "B"
    return condition

run_mapper = lambda x: get_run(x)
condition_mapper = lambda x: get_condition(x)

df.run_id = df.filename.map(run_mapper)
df.filename = df.filename.map(condition_mapper)
df.to_csv("merged_recolumned.tsv", sep = "\t", index=False)


