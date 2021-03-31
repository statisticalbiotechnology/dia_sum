#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 13:26:25 2021

@author: ptruong
"""
import pandas as pd
import numpy as np

import os 

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/")


def prettyPrint(df):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)


df_all = pd.DataFrame()
for file in os.listdir():
    if file[-4:] == ".tsv":
        df = pd.read_csv(file, sep = "\t")
        print(file)
        print(len(df))
        df_all = pd.concat([df_all, df], axis = 0)
df = df_all.reset_index().drop("index", axis = 1)

df = df[df.decoy != 1]
df = df[df.m_score < 0.01]

# filename has different formatting, we need to change number or implement regex.
experiment_id_mapper = lambda x: x.split("_")[5]
sample_id_mapper = lambda x: x.split("_")[8] #hye124 
specie_mapper = lambda x: x.split("_")[-1]

df["experiment_id"] = df["filename"].map(experiment_id_mapper)
df["sample_id"] = df["filename"].map(sample_id_mapper)
df = df[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]
df["specie"] = df.ProteinName.map(specie_mapper)


df.describe()

df.Intensity.min()
df.Intensity.max()

np.log2(df.Intensity).hist(bins = 1000)

experiments_list = df.experiment_id.unique()
specie_list = df.specie.unique()

specie = specie_list[0]
exp1 = experiments_list[0]
exp2 = experiments_list[1]

df_exp1 = df[(df.experiment_id == exp1) & (df.specie == specie)]
df_exp2 = df[(df.experiment_id == exp2) & (df.specie == specie)]


def get_protein_intensity_df(df_exp):
    protein_array = []
    intensity_array = []
    for protein in df_exp.ProteinName.unique():
        intensity = df_exp[df_exp.ProteinName == protein].sort_values(by = "m_score", ascending = True).Intensity.iloc[0]
        protein_array.append(protein)
        intensity_array.append(float(intensity))
        
    df_prot = pd.DataFrame(np.array([protein_array, intensity_array]).T, columns = ["protein", "intensity"])
    if len(df_exp.experiment_id.unique()) == 1:
        df_prot["run"] = df_exp.experiment_id.unique()[0]
    else:
        print("Non-unique run id, this should not happen!")
    if len(df_exp.sample_id.unique()) == 1:
        df_prot["sample"] = df_exp.sample_id.unique()[0]
    else:
        print("Non-unique sample id, this should not happen!")
    df_prot["intensity"] = df_prot["intensity"].astype(float)
    return df_prot    

df_prot1 = get_protein_intensity_df(df_exp1)
df_prot2 = get_protein_intensity_df(df_exp2)


def get_intersecting_protein(df_prot1, df_prot2):
    intersecting_proteins = pd.Series(list(set(df_prot1.protein).intersection(set(df_prot2.protein))))
    return len(intersecting_proteins)
    
def get_differentially_expressed_protein_ECOLI(df_prot1, df_prot2, treshold = 2.0):
    intersecting_proteins = pd.Series(list(set(df_prot1.protein).intersection(set(df_prot2.protein))))
    print(str(len(intersecting_proteins)) + " intersecting proteins.")
    
    fc_array = []
    protein_array = []
    for protein in intersecting_proteins:
        protA = np.log2(df_prot1[df_prot1.protein == protein].intensity).iloc[0]
        protB = np.log2(df_prot2[df_prot2.protein == protein].intensity).iloc[0]
        fc = protB - protA
        fc_array.append(fc)
        protein_array.append(protein)
    df_fc = pd.DataFrame(np.array([protein_array, fc_array]).T, columns = ["protein", "log2_fc"])
    df_fc["log2_fc"] = df_fc["log2_fc"].astype(float)
    return (df_fc["log2_fc"] > treshold).sum()

def get_differentially_expressed_protein_YEAST(df_prot1, df_prot2, treshold = -1.0):
    intersecting_proteins = pd.Series(list(set(df_prot1.protein).intersection(set(df_prot2.protein))))
    print(str(len(intersecting_proteins)) + " intersecting proteins.")
    
    fc_array = []
    protein_array = []
    for protein in intersecting_proteins:
        protA = np.log2(df_prot1[df_prot1.protein == protein].intensity).iloc[0]
        protB = np.log2(df_prot2[df_prot2.protein == protein].intensity).iloc[0]
        fc = protB - protA
        fc_array.append(fc)
        protein_array.append(protein)
    df_fc = pd.DataFrame(np.array([protein_array, fc_array]).T, columns = ["protein", "log2_fc"])
    df_fc["log2_fc"] = df_fc["log2_fc"].astype(float)
    return (df_fc["log2_fc"] < treshold).sum()


def get_differentially_expressed_protein_HUMAN(df_prot1, df_prot2, treshold = 0.5):
    intersecting_proteins = pd.Series(list(set(df_prot1.protein).intersection(set(df_prot2.protein))))
    print(str(len(intersecting_proteins)) + " intersecting proteins.")
    
    fc_array = []
    protein_array = []
    for protein in intersecting_proteins:
        protA = np.log2(df_prot1[df_prot1.protein == protein].intensity).iloc[0]
        protB = np.log2(df_prot2[df_prot2.protein == protein].intensity).iloc[0]
        fc = protB - protA
        fc_array.append(fc)
        protein_array.append(protein)
    df_fc = pd.DataFrame(np.array([protein_array, fc_array]).T, columns = ["protein", "log2_fc"])
    df_fc["log2_fc"] = df_fc["log2_fc"].astype(float)
    return (((df_fc["log2_fc"] > treshold).sum()) + ((df_fc["log2_fc"] < -treshold).sum()))

specie_list

import time
treshold = 2.0
specie = "ECOLI"
start = time.time()
de_array = []
intersecting_array = []
for exp1 in experiments_list:
    de_column_array = []
    intersecting_column_array = []
    print("exp1: " + exp1)
    for exp2 in experiments_list:
        
        df_exp1 = df[(df.experiment_id == exp1) & (df.specie == specie)]
        df_exp2 = df[(df.experiment_id == exp2) & (df.specie == specie)]
        
        df_prot1 = get_protein_intensity_df(df_exp1)
        df_prot2 = get_protein_intensity_df(df_exp2)
        
        intersecting_proteins = get_intersecting_protein(df_prot1, df_prot2)
        if specie == "ECOLI":
            de = get_differentially_expressed_protein_ECOLI(df_prot1, df_prot2, treshold = treshold)
        elif specie == "YEAS8":
            de = get_differentially_expressed_protein_YEAST(df_prot1, df_prot2, treshold = treshold)
        elif specie == "HUMAN":
            de = get_differentially_expressed_protein_HUMAN(df_prot1, df_prot2, treshold = treshold)
 
        intersecting_column_array.append(intersecting_proteins)
        de_column_array.append(de)
        print(time.time() - start)
        print("exp2: " + exp2)

    de_array.append(de_column_array)
    intersecting_array.append(intersecting_column_array)   
    
exp_group_list = []
for exp in experiments_list:
    group = df[df.experiment_id == exp].sample_id.unique()[0] 
    exp_sample = exp + "_group_" + group
    exp_group_list.append(exp_sample)
    
end = time.time()
print(end-start)

df_de = pd.DataFrame(de_array, columns = [exp + "_log2(B)" for exp in exp_group_list], 
                     index = [exp + "_log2(A)" for exp in exp_group_list])


prettyPrint(df_de)

from IPython.display import display

display(df_de.style.bar(color='#d65f5f'))


                        
############
# Triqler ## 
############
os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_output_to_df import parse_triqler

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_ecoli_fc2_full_tsv")

df = parse_triqler("proteins.tsv")
df = df[df.protein_id_posterior_error_prob < 0.05]
specie_mapper = lambda x : x.split("_")[-1]
df["specie"] = df.protein.map(specie_mapper)

specie = "ECOLI"
treshold = 2.0

    
def get_triqler_DE_fc(df, specie, treshold):
    df_rv = df[df.specie == specie].drop(['protein', 'peptides', 'q_value', 'posterior_error_prob',
           'num_peptides', 'protein_id_posterior_error_prob', 'log2_fold_change',
           'diff_exp_prob_2.0', "specie"], axis = 1 ).astype(float)
        
    def get_fc_DE_ECOLI(A, B, treshold = 2.0):
        return (np.log2(B / A) > treshold).sum()
    def get_fc_DE_YEAST(A, B, treshold = -1.0):
        return (np.log2(B / A) < treshold).sum()
    def get_fc_DE_HUMAN(A, B, treshold = 0.5):
        return ((np.log2(B / A) > treshold).sum() + 
                (np.log2(B / A) < -treshold).sum())
    
    de_array = []
    start = time.time()
    for A in df_rv.columns:
        de_column_array = []
        print("A: " + A)
    
        for B in df_rv.columns:
            print("B: " + B)
            if specie == "HUMAN":
                de = get_fc_DE_HUMAN(df_rv[A],df_rv[B], treshold = treshold)
            elif specie == "YEAS8":
                de = get_fc_DE_YEAST(df_rv[A],df_rv[B], treshold = treshold)
            elif specie == "ECOLI":
                de = get_fc_DE_ECOLI(df_rv[A],df_rv[B], treshold = treshold)
            de_column_array.append(de)
        de_array.append(de_column_array)
        print(time.time() - start)
            
    
    end = time.time()
    print(end-start)
    
    df_de = pd.DataFrame(de_array, columns = [exp + "_log2(B)" for exp in df_rv.columns], 
                         index = [exp + "_log2(A)" for exp in df_rv.columns])
    return df_de
        

A = "1:002-Pedro"
B = "2:003-Pedro"


