#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:35:02 2021

@author: ptruong
"""

#import pandas as pd 
#import numpy as np 

def fasta_to_protein_specie_map(filename):
    f = open(filename, "r")
    #line = f.readline() 
    
    proteins = []
    species = []
    for line in f:
        if line[0] == ">":
            protein = line.split("|")[1]
            specie = line.split("|")[2].split(" ")[0].split("_")[1]
            proteins.append(protein)
            species.append(specie)
    
    #protein_specie_map = pd.DataFrame([proteins, species], index = ["protein", "specie"]).T
    return dict(zip(proteins, species))

if __name__ == "__main__":
    filename = "2021-04-27-decoys-reviewed-contam-UP000005640-UP000002311-UP000000625.fas"
    #protein_specie_map.to_csv("protein_specie_map.csv", sep = "\t", index = False)
    # test to see if species have same protein... they don't seem to have it.
    #for protein in protein_specie_map.protein:
    #    if len(protein_specie_map[protein_specie_map.protein==protein].specie.unique()) > 1:
    #        print(protein_specie_map[protein_specie_map.protein==protein])

