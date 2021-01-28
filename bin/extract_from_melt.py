#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 00:55:45 2020

@author: ptruong
"""

import pandas as pd
import numpy as np

def get_protein_abundance(df, sample, specie, protein):
    """
    Get the protein abundances from all runs for a sample and specie.
    
    df = triq or spec (melted format)
    sample = S01...S10
    specie = {HUMAN; ARATH; CAEEL}
    """
    df_sample = df[df["sample"] == sample]
    df_sample_specie = df_sample[df_sample["specie"] == specie]
    df_sample_specie_protein = df_sample_specie[df_sample_specie["protein"] == protein]
    
    return df_sample_specie_protein.value

def get_proteins(df, specie):
    df_specie = df[df["specie"] == specie]
    #df_specie_protein = df_specie[df_specie["protein"] == protein]
    proteins = df_specie.protein.unique()
    return proteins    











