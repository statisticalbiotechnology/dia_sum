#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 12:03:09 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

from constants import get_log2FC_ratio_matrices

def get_differentially_expressed_proteins_from_log2FC_df(method, sample1, sample2, specie, ratio, two_sided = False):
    """
    Takes log2FC df from get_log2FC_matrix() and computes differentially expressed 
    proteins for each run.
    
    example params:
    
        method = "triq" # VAR
        sample1 = "S02" # VAR
        sample2 = "S06" # VAR
        specie = "HUMAN" # VAR
        ratio = 0.8 # VAR - ratio for what is considered diff. exp.
        two_sided = False # VAR - is this needed? I will skip this for now
    
    NOTE: two_sided does nothing atm. left for furutre fix.
    """
    ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix = get_log2FC_ratio_matrices()
    
    if specie == "HUMAN":
        FC_matrix = HUMAN_FC_matrix
        FC_treshold = ratio * FC_matrix[sample2][sample1] #The smaples are reversed because the log2FC in this function is reversed.
        file_name = method+"_"+"log2FC"+"_"+sample1+"_"+sample2+"_"+specie+".csv"
        df_log2FC = pd.read_csv(file_name, sep = "\t")
        diffExp = np.sum(df_log2FC[df_log2FC > FC_treshold] > 0) #This comparison should be more then
    elif specie == "ARATH":
        FC_matrix = ARATH_FC_matrix
        FC_treshold = ratio * FC_matrix[sample2][sample1] #The smaples are reversed because the log2FC in this function is reversed.
        file_name = method+"_"+"log2FC"+"_"+sample1+"_"+sample2+"_"+specie+".csv"
        df_log2FC = pd.read_csv(file_name, sep = "\t")
        diffExp = np.sum(df_log2FC[df_log2FC < -ratio] > 0) + np.sum(df_log2FC[df_log2FC > ratio] > 0) #This comparison should be if less or more than ratio
    elif specie == "CAEEL":
        FC_matrix = CAEEL_FC_matrix
        FC_treshold = ratio * FC_matrix[sample2][sample1] #The smaples are reversed because the log2FC in this function is reversed.
        file_name = method+"_"+"log2FC"+"_"+sample1+"_"+sample2+"_"+specie+".csv"
        df_log2FC = pd.read_csv(file_name, sep = "\t")
        diffExp = np.sum(df_log2FC[df_log2FC < FC_treshold] > 0) #This comparison should be less than
    else:
        print("no species specificed")
        # return 
    return diffExp
