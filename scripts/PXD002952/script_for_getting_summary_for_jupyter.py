#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 23:42:24 2021

@author: ptruong
"""

# A simple script for the notebook just to simply


import itertools


import pandas as pd
import numpy as np


def get_summary_hye110():
    df = pd.read_csv("aligned.csv", sep = "\t")
    
    experiment_id_mapper = lambda x: x.split("_")[5]
    sample_id_mapper = lambda x: x.split("_")[7]
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
    
    return df_summary    
    
    
def get_summary_hye124():
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
    
    return df_summary    
    
    
   
def get_summary_hye124_ttof6600_64var():
    df = pd.read_csv("aligned.csv", sep = "\t")
    
    experiment_id_mapper = lambda x: x.split("_")[5]
    sample_id_mapper = lambda x: x.split("_")[8]
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
    
    return df_summary    
    

#################
# Get diff exp ##
#################
