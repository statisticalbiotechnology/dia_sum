#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:21:30 2020

@author: ptruong
"""

import pandas as pd
import numpy as np 

def normalize_within_sample(df):
    """
    After this every sample and run will have protein abundances sum to 1.
    i.e. S01 R01 will have proteins for CAEEL,HUMAN and ARATH sum to 1. This
    should give us thier relative abundances within smaple.
    
    df is the melted df from 
        
    triq = melt_triqler_output(triq)
    spec = melt_spectronaut_triqler_formatted(spec)
    
    """
    samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 
    runs = ["R0" + str(i) for i in range(1,6)]
    
    df_norm = pd.DataFrame()
    for i in samples:
        for j in runs:
            df_sample = df[df["sample"] == i]
            df_sample_run = df_sample[df_sample["run"] == j]
            
            df_sample_run.value = df_sample_run.value / df_sample_run.value.sum() # Normalize
            df_norm = df_norm.append(df_sample_run)
    return df_norm

def get_ratios_from_normalized_melted_df(df, specie):
    """
    input - normalized df (triq, spec)
    output - ratios dataframe
    """
    samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 
    runs = ["R0" + str(i) for i in range(1,6)]
    
    ratios = []
    for sample in samples:
        ratios_for_sample = []
        for run in runs:
            df_sample = df[df["sample"] == sample]
            df_sample_run = df_sample[df_sample["run"] == run]
            df_sample_run_spec = df_sample_run[df_sample_run["specie"] == specie]
            ratios_for_sample.append(df_sample_run_spec.value.sum())
        ratios.append(ratios_for_sample)
    ratios_df = pd.DataFrame(ratios, index = samples, columns = runs)
    return ratios_df


