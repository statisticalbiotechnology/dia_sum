#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:15:19 2021

@author: ptruong
"""


##########
# OSW DE #
##########

import pandas as pd
import numpy as np

import os 

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/full_ts_v")


def prettyPrint(df):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)

# Readin and filter 
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

# Check intensity histogram for batch effect and noise tail.
np.log2(df.Intensity).hist(bins = 1000)

experiments_list = df.experiment_id.unique()
specie_list = df.specie.unique()


# Select the protein intensity with highest searchScore, i.e. lowest m_score, same as triqler does, but we put
# in -np.log10(m_score) as triqler search results.
def get_protein_intensity_df(df_exp):
    protein_array = []
    intensity_array = []
    for protein in df_exp.ProteinName.unique():
        intensity = df_exp[df_exp.ProteinName == protein].sort_values(by = "m_score", ascending = False).Intensity.iloc[0]
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

# Defining some functions for finding differential expressions for different species.
# Defining function for finding intersecting proteins between two samples.

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



import time
#treshold = 2.0
#specie = "ECOLI"

def get_DE_df(df, specie, treshold):
    start = time.time()
    de_array = []
    intersecting_array = []
    for exp1 in experiments_list:
        de_column_array = []
        intersecting_column_array = []
        print("")
        print("A: " + exp1)
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
            print(time.time() - start, end = "")
            print("s")
            print("B: " + exp2, end = " : ")

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
    return df_de


def get_mean_intensities_per_sample(df):
    samples = df.sample_id.unique()
    dfs = {}
    for sample in samples:
        proteins = df[df.sample_id == sample].ProteinName.unique()
        df_sample = df[df.sample_id == sample]
        
        mu_array = []
        std_array= []
        for protein in proteins:
            mu_intensity = df_sample[df_sample["ProteinName"] == protein].Intensity.mean()
            std_intensity = df_sample[df_sample["ProteinName"] == protein].Intensity.std()
            mu_array.append(mu_intensity)
            std_array.append(std_intensity)
        df_intensity = pd.DataFrame(np.array([mu_array, std_array]).T, index = proteins, columns = ["mu", "std"])
        dfs[sample] = df_intensity
    return dfs

def get_intersecting_protein(df_prot1, df_prot2, protein_index = False):
    if protein_index == True:
        intersecting_proteins = pd.Series(list(set(df_prot1.index).intersection(set(df_prot2.index))))
    else:
        intersecting_proteins = pd.Series(list(set(df_prot1.protein).intersection(set(df_prot2.protein))))
    return len(intersecting_proteins)

dfs = get_mean_intensities_per_sample(df)
get_intersecting_protein(dfs["1"], dfs["2"], protein_index = True)



dfs["1"]






df_de_ECOLI = get_DE_df(df, specie = "ECOLI", treshold = 2)

specie_list
dfs = []


for i in range(0,200+5,5):
    i = i/100
    de = get_DE_df(df, specie = "ECOLI", treshold = i)
    
##############
# Triqler DE #
##############

import os 

import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

from triqler_output_to_df import  parse_triqler
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_results")


de_array = []
for file in sorted(os.listdir()):
    df = parse_triqler(file)
    de_protein = len(df[df.q_value < 0.05])
    de_array.append(de_protein)




