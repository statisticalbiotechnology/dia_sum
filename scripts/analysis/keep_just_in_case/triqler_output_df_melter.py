#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 17:09:04 2020

@author: ptruong
"""

import numpy as np
import pandas as pd

import time


# def function for melt spec

def melt_spectronaut_triqler_formatted(spec):
    """
    Input
        spec - triqler formatted spectronaut data. (same as .pkl file)
    """
    spec = spec.rename(columns={'S03:S04_R05': 'S03:S03_R04'})
    sample_run_cols  = ['S01:S01_R01', 'S01:S01_R02', 'S01:S01_R03',
           'S01:S01_R04', 'S01:S01_R05', 'S02:S02_R01', 'S02:S02_R02',
           'S02:S02_R03', 'S02:S02_R04', 'S02:S02_R05', 'S03:S03_R01',
           'S03:S03_R02', 'S03:S03_R03', 'S03:S03_R05', 'S03:S03_R04',
           'S04:S04_R01', 'S04:S04_R02', 'S04:S04_R03', 'S04:S04_R04',
           'S04:S04_R05', 'S05:S05_R01', 'S05:S05_R02', 'S05:S05_R03',
           'S05:S05_R04', 'S05:S05_R05', 'S06:S06_R01', 'S06:S06_R02',
           'S06:S06_R03', 'S06:S06_R04', 'S06:S06_R05', 'S07:S07_R01',
           'S07:S07_R02', 'S07:S07_R03', 'S07:S07_R04', 'S07:S07_R05',
           'S08:S08_R01', 'S08:S08_R02', 'S08:S08_R03', 'S08:S08_R04',
           'S08:S08_R05', 'S09:S09_R01', 'S09:S09_R02', 'S09:S09_R03',
           'S09:S09_R04', 'S09:S09_R05', 'S10:S10_R01', 'S10:S10_R02',
           'S10:S10_R03', 'S10:S10_R04', 'S10:S10_R05']
    
    data = spec
    cols = list(data.drop(sample_run_cols, axis = 1).columns) + ["sample", "run", "value"]
    start = time.time()
    vals = []
    for i in range(len(data)):
        #print(str(i)+ "/" + str(len(data)))    
        id_val = data.iloc[i,:][0]
        specie_val = data.iloc[i,:][1]
        protein_val = data.iloc[i,:][2]
        for j in range(len(data.iloc[i,:])):
            #loop through data, and labels to melt the df#
            col =  data.iloc[i,:].index[j]
            if col[0] == "S":
                sample, run = col.split(":")[1].split("_")
                val = data.iloc[i,:][j]
                vals.append([id_val, specie_val, protein_val, sample, run, val])
    end = time.time()
    #print(end-start)
    
    melted = pd.DataFrame(vals, columns = cols)
    return melted



def melt_triqler_output(triq, protein_id_pep_treshold = 0.01):
    # treshold PEP triq
    triq = triq[triq.protein_id_posterior_error_prob < protein_id_pep_treshold]
    
    # def cuntion for melt triq
    sample_run_cols = ['S01:S01_R01', 'S01:S01_R02', 'S01:S01_R03',
           'S01:S01_R04', 'S01:S01_R05', 'S02:S02_R01', 'S02:S02_R02',
           'S02:S02_R03', 'S02:S02_R04', 'S02:S02_R05', 'S03:S03_R01',
           'S03:S03_R02', 'S03:S03_R03', 'S03:S03_R04', 'S03:S03_R05',
           'S04:S04_R01', 'S04:S04_R02', 'S04:S04_R03', 'S04:S04_R04',
           'S04:S04_R05', 'S05:S05_R01', 'S05:S05_R02', 'S05:S05_R03',
           'S05:S05_R04', 'S05:S05_R05', 'S06:S06_R01', 'S06:S06_R02',
           'S06:S06_R03', 'S06:S06_R04', 'S06:S06_R05', 'S07:S07_R01',
           'S07:S07_R02', 'S07:S07_R03', 'S07:S07_R04', 'S07:S07_R05',
           'S08:S08_R01', 'S08:S08_R02', 'S08:S08_R03', 'S08:S08_R04',
           'S08:S08_R05', 'S09:S09_R01', 'S09:S09_R02', 'S09:S09_R03',
           'S09:S09_R04', 'S09:S09_R05', 'S10:S10_R01', 'S10:S10_R02',
           'S10:S10_R03', 'S10:S10_R04', 'S10:S10_R05']
    data = triq
    cols = list(data.drop((sample_run_cols + ["peptides"]), axis = 1).columns) + ["sample", "run", "value", "peptide"]
    
    
    
    start = time.time()
    vals = []
    for i in range(len(data)):
       #print(str(i)+ "/" + str(len(data)))    
        q_val = data.iloc[i,:][0]
        posterior_error_prob_val = data.iloc[i,:][1]
        protein_val = data.iloc[i,:][2]
        num_peptides_val = data.iloc[i,:][3]
        protein_id_posterior_error_prob_val = data.iloc[i,:][4]
        log2_fold_change_val = data.iloc[i,:][5]
        diff_exp_prob_val = data.iloc[i,:][6]
        peptides_val = data.iloc[i,:][-1]
        for j in range(len(data.iloc[i,:])):
            #loop through data, and labels to melt the df#
            col =  data.iloc[i,:].index[j]
            if col[0] == "S":
                sample, run = col.split(":")[1].split("_")
                val = data.iloc[i,:][j]
                vals.append([q_val, posterior_error_prob_val, protein_val, num_peptides_val,
                             protein_id_posterior_error_prob_val, log2_fold_change_val,
                             diff_exp_prob_val, sample, run, val, peptides_val])
    end = time.time()
    #print(end-start)
    
    melted = pd.DataFrame(vals, columns = cols)
    melted["specie"] = melted.protein.apply(lambda x: x.split('_')[-1])
    return melted
    
    
    
