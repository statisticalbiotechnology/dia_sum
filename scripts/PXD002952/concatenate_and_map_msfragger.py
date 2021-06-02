#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 10:35:30 2021

@author: ptruong
"""

import os

os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger") 

import pandas as pd
import re


sample_mapper = re.findall(file)

def sample_mapper(file): return re.findall(r"(?<=Sample_)\d+", file)[0] 
def replicate_mapper(file): return re.findall(r"(?<=Repl)\d+", file)[0]
def quality_level_mapper(file): return re.findall(r"(?<=_Q)\d+", file)[0]
def run_mapper(file): return re.findall(r"\d+(?=-Pedro)", file)[0]

def merge_msfragger():
    df = pd.DataFrame()
    for file in os.listdir():
        if file[:8] == "interact":
            if file[-4:] == ".tsv":
                print(file)
                try:
                    df_ = pd.read_csv(file, sep = "\t")
                    df_["sample"] = df_.spectrum.map(sample_mapper)
                    df_["replicate"] = df_.spectrum.map(replicate_mapper)
                    df_["quality_level"] = df_.spectrum.map(quality_level_mapper)
                    df_["run"] = df_.spectrum.map(run_mapper)
                    df = pd.concat([df, df_], axis = 0)
                except:
                    print("FAIL: " + file)
    return df

df = merge_msfragger()
# Select only top quality spectras
df = df[df["quality_level"] == "1"]

# Write triqler converter
df_triqler = pd.DataFrame()
df_triqler["run"] = df["run"]
df_triqler["condition"] = df["sample"]
df_triqler["charge"] = df["assumed_charge"]
df_triqler["searchScore"] = df["peptideprophet_probability"]
df_triqler["intensity"] = df[""]
df_triqler["peptide"] = df[""]
df_triqler["proteins´"] =  df[""]




