#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 23:18:34 2021

@author: ptruong
"""
import itertools
import re

import pandas as pd
import numpy as np

df = pd.read_csv("aligned.csv", sep = "\t")

experiment_id_mapper = lambda x: x.split("_")[5]
sample_id_mapper = lambda x: x.split("_")[9]
df["experiment_id"] = df["filename"].map(experiment_id_mapper)
df["sample_id"] = df["filename"].map(sample_id_mapper)

df = df[df.m_score < 0.01] #FDR-filtering
df = df[df.decoy == 0] #remove decoy

def get_protein_quantity_df(df_sample):
    """
    E.g. input of df_sample:
        
    df_sample = df[df.experiment_id == "001-Pedro"]
    """
    def get_protein_quantity(protein):
        df_prot = df_sample[df_sample.ProteinName == protein]
        if len(df_prot) > 2:
            intensity = df_prot.Intensity.sort_values(ascending=False).head(3).mean()
        else:
            intensity = np.nan
        return intensity
    
    import time
    start=time.time()
    df_sample["proteinQuantity"] = df_sample["ProteinName"].map(get_protein_quantity)
    end=time.time()
    print("time: " + str(end-start) + " s")
    prot = df_sample[df_sample["proteinQuantity"].notna()]
    prot = prot[["experiment_id", "sample_id", "ProteinName", "proteinQuantity"]]
    #prot = prot.drop_duplicates()
    return prot


def protein_quantity_df(df):
    """
    df is aligned.csv from TRIC, with experiment_id and smaple_id mapping. 
    
    """
    df_protein = pd.DataFrame()
    for sample_id in df.experiment_id.unique():
        df_sample = df[df.experiment_id == sample_id]    
        df_protein_subset = get_protein_quantity_df(df_sample)
        df_protein = df_protein.append(df_protein_subset)
    return df_protein


def get_summary(df_protein):
    """
    df_protein is the output from protein_quantity_df
    """
    exp_array = []
    sample_array = []
    proteins_array = []
    for exp in df_protein.experiment_id.unique():
        df_exp = df_protein[df_protein.experiment_id == exp]
        sample = df_exp.sample_id.unique()
        n_proteins = len(df_exp)
        exp_array.append(exp)
        sample_array.append(sample)
        proteins_array.append(n_proteins)
    df_summary = pd.DataFrame(np.array([exp_array, sample_array, proteins_array]).T, columns = ["experiment_id", "sample_id", "proteins"])
    return df_summary

df_prot = protein_quantity_df(df) # Protein quantification, duplicates_dropped.

df_summary = get_summary(df_prot) # Check how many proteins quantified in each sample.
df_summary["peptides"] = df_summary.proteins 
df_summary["proteins"] = get_summary(df_prot.drop_duplicates()).proteins

df_prot = df_prot.drop_duplicates()

# Check overlap of protein between samples (?needed)
# Check diff expression of overlapping proteins.


def get_sample_by_sample_proteinQuantities(df_prot):
    def clean_df_x(df_x):
        """
        cleanup for experiment df
        """
        experiment_id_x = df_x.experiment_id.unique()[0]
        sample_id_x = df_x.sample_id.unique()[0]
        df_x = df_x.drop(["experiment_id", "sample_id"], axis = 1)
        df_x = df_x.rename({"proteinQuantity": experiment_id_x + "_" + sample_id_x + "_proteinQuantity"}, axis = 1)
    
        return df_x
    
    df_x_list = []
    for exp in df_prot.experiment_id.unique():
        df_x = df_prot[df_prot.experiment_id == exp].set_index("ProteinName")
        df_x = clean_df_x(df_x)
        df_x_list.append(df_x)
    df_xs = pd.concat(df_x_list, axis = 1)
    
    return df_xs

def proteinName_to_specie_mapper(x):
    specie_set = ['YEAS8', 'HUMAN', 'ECOLI']
    protein_str_as_list = re.split("/|_", "_".join(x.split("|")))
    
    specie_list = list(set(specie_set).intersection(set(protein_str_as_list)))
    if len(specie_list) > 1:
        specie = x
    else:
        specie = specie_list[0]
    return specie


def get_all_comps(sample_set):
    """
    Get comparisons of samples lists.
    """
    comps = []
    for result in itertools.combinations(sample_set, 2):
        comps.append(result)
    return comps


df_xs = get_sample_by_sample_proteinQuantities(df_prot)
sample_set = list(df_xs.columns)
comps = get_all_comps(sample_set)
df_xs["specie"] = df_xs.index.map(proteinName_to_specie_mapper)

#### compute differential expression - sample by sample



def get_tresholded_HUMAN_diffExp(df_xs, exp_1, exp_2, treshold = 0.5):
    """
    df_xs is output from get_sample_by_sample_proteinQuantities(df_prot) 
    """
    df_1 = df_xs[exp_1]
    df_2 = df_xs[exp_2]
    
    df_diff = pd.DataFrame(np.log2(df_1)-np.log2(df_2), columns = ["log2_proteinQuantity"])
    df_diff = df_diff.dropna()
    df_diff["specie"] = df_diff.index.map(proteinName_to_specie_mapper)

    df_diff_spec = pd.DataFrame(df_diff[df_diff.specie == "HUMAN"].log2_proteinQuantity)
    df_treshold_diff_spec = df_diff_spec[(df_diff_spec.log2_proteinQuantity < treshold) & (df_diff_spec.log2_proteinQuantity > - treshold)]
    return df_treshold_diff_spec

def get_tresholded_ECOLI_diffExp(df_xs, exp_1, exp_2, treshold = 2):
    """
    df_xs is output from get_sample_by_sample_proteinQuantities(df_prot) 

    log2(exp_1) - log2(exp_2), so have the larger exp as exp_1...
    """
    
    df_1 = df_xs[exp_1]
    df_2 = df_xs[exp_2]
    
    df_diff = pd.DataFrame(np.log2(df_1)-np.log2(df_2), columns = ["log2_proteinQuantity"])
    df_diff = df_diff.dropna()
    df_diff["specie"] = df_diff.index.map(proteinName_to_specie_mapper)
    
    df_diff_spec = pd.DataFrame(df_diff[df_diff.specie == "ECOLI"].log2_proteinQuantity)
    df_treshold_diff_spec = df_diff_spec[df_diff_spec > treshold]
    return df_treshold_diff_spec

def get_tresholded_YEAST_diffExp(df_xs, exp_1, exp_2, treshold = 2):
    """
    df_xs is output from get_sample_by_sample_proteinQuantities(df_prot) 

    log2(exp_1) - log2(exp_2), so have the larger exp as exp_1...
    """
    
    df_1 = df_xs[exp_1]
    df_2 = df_xs[exp_2]
    
    df_diff = pd.DataFrame(np.log2(df_1)-np.log2(df_2), columns = ["log2_proteinQuantity"])
    df_diff = df_diff.dropna()
    df_diff["specie"] = df_diff.index.map(proteinName_to_specie_mapper)
    
    df_diff_spec = pd.DataFrame(df_diff[df_diff.specie == "YEAS8"].log2_proteinQuantity)
    df_treshold_diff_spec = df_diff_spec[df_diff_spec > treshold]
    return df_treshold_diff_spec


#comp = comps[0]
#human_treshold = 0.2
#get_tresholded_HUMAN_diffExp(df_xs, comp[0], comp[1], treshold = 0.5).dropna()
#get_tresholded_ECOLI_diffExp(df_xs, comp[1], comp[0], treshold = 3.5).dropna()
#get_tresholded_YEAST_diffExp(df_xs, comp[0], comp[1], treshold = 3.5).dropna()

def get_de(df_xs, human_treshold = 0.1, ecoli_treshold = 3.5, yeast_treshold = 1.5):
    def get_HUMAN_de(df_xs, treshold = 0.1):
        comps_array = []
        human_de_array = []
        for comp in comps:
            human_de = len(get_tresholded_HUMAN_diffExp(df_xs, comp[0], comp[1], treshold = treshold).dropna())
            comps_array.append(comp)
            human_de_array.append(human_de)
        return pd.DataFrame(human_de_array, index = comps_array, columns = ["human_de"])
    
    #think about arranging comp
    
    def get_ECOLI_de(df_xs, treshold = 3.5):
        comps_array = []
        ecoli_de_array = []
        for comp in comps:
            if comp[0].split("_")[1] == "A":
                ecoli_de = len(get_tresholded_ECOLI_diffExp(df_xs, comp[1], comp[0], treshold = treshold).dropna())
                comps_array.append(comp)
                ecoli_de_array.append(ecoli_de)
            else:
                ecoli_de = len(get_tresholded_ECOLI_diffExp(df_xs, comp[0], comp[1], treshold = treshold).dropna())
                comps_array.append(comp)
                ecoli_de_array.append(ecoli_de)       
        return pd.DataFrame(ecoli_de_array, index = comps_array, columns = ["ecoli_de"])
    
    def get_YEAST_de(df_xs, treshold = 1.5):
        comps_array = []
        yeast_de_array = []
        for comp in comps:
            if comp[0].split("_")[1] == "A":
                yeast_de = len(get_tresholded_YEAST_diffExp(df_xs, comp[0], comp[1], treshold = treshold).dropna())
                comps_array.append(comp)
                yeast_de_array.append(yeast_de)
            else:
                yeast_de = len(get_tresholded_YEAST_diffExp(df_xs, comp[1], comp[0], treshold = treshold).dropna())
                comps_array.append(comp)
                yeast_de_array.append(yeast_de)       
        return pd.DataFrame(yeast_de_array, index = comps_array, columns = ["yeast_de"])
    
    human_de = get_HUMAN_de(df_xs, treshold = human_treshold)
    ecoli_de = get_ECOLI_de(df_xs, treshold = ecoli_treshold)
    yeast_de = get_YEAST_de(df_xs, treshold = yeast_treshold)
    return pd.concat([human_de, ecoli_de, yeast_de], axis = 1)

de = get_de(df_xs, human_treshold = 0.1, ecoli_treshold = 3.5, yeast_treshold = 1.5)


