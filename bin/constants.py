#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 00:39:28 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

def get_sample_ratios():
    ARATH = np.array([0.5,  0.5,  0.5,  0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    CAEEL = np.array([ 0.5,  0.25,  0.125,  0.0625,  0.031,  0.0155, 0.008,  0.004,  0.002,  0.000001])
    HUMAN = np.array([ 00.000001,  0.25, 0.375,  0.4375,  0.469,  0.4845,  0.492, 0.496, 0.498,  0.5])
    return ARATH, CAEEL, HUMAN

def get_log2FC_ratio_matrix(specie_array):
    """
    Inpur specie array with mixture ratios, and generate log2FC matrix.
    """
    samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 
    FC_ratios = []
    for i in range(10):
        sample_1 = specie_array[i]
        FC_ratio = []
        for j in range(10):
            sample_2 = specie_array[j]
            FC_ratio.append(np.log2(sample_2) - np.log2(sample_1))
        FC_ratios.append(FC_ratio)
    df_FC_ratios = pd.DataFrame(FC_ratios, index = samples, columns = samples)
    return df_FC_ratios

def get_log2FC_ratio_matrices():
    ARATH = np.array([0.5,  0.5,  0.5,  0.5,  0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    CAEEL = np.array([ 0.5,  0.25,  0.125,  0.0625,  0.031,  0.0155, 0.008,  0.004,  0.002,  0.000001])
    HUMAN = np.array([ 00.000001,  0.25, 0.375,  0.4375,  0.469,  0.4845,  0.492, 0.496, 0.498,  0.5])
    
    ARATH_FC_matrix = get_log2FC_ratio_matrix(ARATH)
    CAEEL_FC_matrix = get_log2FC_ratio_matrix(CAEEL)
    HUMAN_FC_matrix = get_log2FC_ratio_matrix(HUMAN)
    return ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix
