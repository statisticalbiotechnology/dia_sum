#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 08:16:17 2022

@author: ptruong
"""

import pandas as pd 
import numpy as np

#df = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_002-Pedro_-_Sample_1_-_SW32_-_Repl1.mzML_with_dscore.csv", sep = "\t")
df = pd.read_csv("concatenated_osw_results.csv", sep = "\t")


def map_zero_transitions(value, threshold):
    """
    Use aggr_Peak_Area to map zero/low intensity transition. 
    """
    x = value.split(";")
    x = np.array(x).astype(float)
    index = np.where(x < threshold)[0]    
    mask = np.ones(x.shape[0], dtype=bool)
    mask[index] = False    
    return mask

def map_top_n_transitions(x, n):    
    x = x.split(";")
    x = np.array(x).astype(float)
    top_n_indices = x.argsort()[-n:][::-1]
    mask = np.zeros(x.shape[0], dtype=bool)
    mask[top_n_indices] = True 
    return mask

def mask_transitions(x, mask):
    x = np.array(x.split(";"))
    return x[mask]


def filter_transition(df, filter_function, *args):
    """
    Function to keep only transition over treshold intensity.
    """
    df["transition_mask"] = df["aggr_Peak_Area"].map(lambda x: filter_function(x, args[0]))
    
    df["aggr_Peak_Area"] = df.apply(lambda x: mask_transitions(x["aggr_Peak_Area"], x["transition_mask"]), axis = 1)
    df["aggr_Peak_Apex"] = df.apply(lambda x: mask_transitions(x["aggr_Peak_Apex"], x["transition_mask"]), axis = 1)
    df["aggr_Fragment_Annotation"] = df.apply(lambda x: mask_transitions(x["aggr_Fragment_Annotation"], x["transition_mask"]), axis = 1)
    df.drop("transition_mask", axis = 1, inplace = True)
    
    df["aggr_Peak_Area"] = df["aggr_Peak_Area"].map(lambda x: ";".join(x))
    df["aggr_Peak_Apex"] = df["aggr_Peak_Apex"].map(lambda x: ";".join(x))
    df["aggr_Fragment_Annotation"] = df["aggr_Fragment_Annotation"].map(lambda x: ";".join(x))

    return df

#df_ = filter_transition(df, map_zero_transitions, 0.05)  # removes low intensity transitions.
df_ = filter_transition(df, map_top_n_transitions, 6) # selects top 6 transitions. (Similar to DIA-NN)
df_.drop("Unnamed: 0", axis = 1, inplace = True)
df_.to_csv("concatenated_osw_results_transitions_filtered_n_6.csv", sep = "\t", index = False)


