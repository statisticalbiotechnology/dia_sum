#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 01:21:59 2020

@author: ptruong
"""

os.chdir("/home/ptruong/git/bayesMS/bin")

import os 

import pandas as pd
import numpy as np 

from read_triqler_output import read_triqler_protein_output_to_df
from triqler_output_df_melter import melt_spectronaut_triqler_formatted, melt_triqler_output


spec = pd.read_pickle(r'spectronaut.pkl')
spec = spec.rename(columns={'S03:S04_R05': 'S03:S03_R04'})
triq = pd.read_pickle(r'triqler.pkl')

triq = read_triqler_protein_output_to_df('proteins.3vs8.tsv')


triq = melt_triqler_output(triq)
spec = melt_spectronaut_triqler_formatted(spec)


samples = ["S0"+str(i) for i in range(1,10)] + ["S10"]


triq.columns

triq[triq["protein_id_posterior_error_prob"] > 0.01] #Both triqler and spec are already 1% FDR tresholded.


# Identified protein through all samples
print("%s : %i" % ("triqler number of protein ids", int(sum(triq.protein_id_posterior_error_prob < 0.01))))
print("%s : %i" % ("spectronaut number of protein ids", int(len(spec.dropna()))))

# Identified protein for species across all samples
triq_count = len(triq[triq.specie == "ARATH"])
spec_count = len(spec[spec.specie == "ARATH"].dropna())

print("%s : %i" % ("triqler ARATH number of protein ids", triq_count))
print("%s : %i" % ("spectronaut ARATH number of protein ids", spec_count))

triq_count = len(triq[triq.specie == "CAEEL"])
spec_count = len(spec[spec.specie == "CAEEL"].dropna())

print("%s : %i" % ("triqler CAEEL number of protein ids", triq_count))
print("%s : %i" % ("spectronaut CAEEL number of protein ids", spec_count))

triq_count = len(triq[triq.specie == "HUMAN"])
spec_count = len(spec[spec.specie == "HUMAN"].dropna())

print("%s : %i" % ("triqler HUMAN number of protein ids", triq_count))
print("%s : %i" % ("spectronaut HUMAN number of protein ids", spec_count))

# Create an identification matrix (df) for each all samples, that can have species as input

samples

specie = None
def count_melted_for_all_samples(df, specie = None):
    """
    df = melted triq or spec.
    """
    counts = []
    for i in samples:
        count_df = df[df["sample"] == i]
        if specie == "ARATH":
            count_df = count_df[count_df["specie"] == "ARATH"].dropna()
        elif specie == "HUMAN":
            count_df = count_df[count_df["specie"] == "HUMAN"].dropna()
        elif specie == "CAEEL":
            count_df = count_df[count_df["specie"] == "CAEEL"].dropna()
        count = len(count_df.dropna())
        counts.append(count)
    return counts

count_melted_for_all_samples(triq, specie = None)
count_melted_for_all_samples(triq, specie = "HUMAN")
count_melted_for_all_samples(triq, specie = "CAEEL")
count_melted_for_all_samples(triq, specie = "ARATH")

count_melted_for_all_samples(spec, specie = None)
count_melted_for_all_samples(spec, specie = "HUMAN")
count_melted_for_all_samples(spec, specie = "CAEEL")
count_melted_for_all_samples(spec, specie = "ARATH")



# Normalization

from normalize_melted import normalize_within_sample, get_ratios_from_normalized_melted_df

df_spec = normalize_within_sample(spec)
df_triq  = normalize_within_sample(triq)

# Explore normalization


get_ratios_from_normalized_melted_df(df_triq, "HUMAN")
get_ratios_from_normalized_melted_df(df_triq, "ARATH")
get_ratios_from_normalized_melted_df(df_triq, "CAEEL")



# Compute True FC on samples

from constants import get_sample_ratios, get_log2FC_ratio_matrix, get_log2FC_ratio_matrices


ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix = get_log2FC_ratio_matrices()
    
# Computing FC diff-exp for triq and spec

from extract_from_melt import get_protein_abundance, get_proteins

proteins = get_proteins(spec, "HUMAN") #FUNC TEST
proteins = get_proteins(triq, "HUMAN") #FUNC TEST
get_protein_abundance(triq, "S01", "HUMAN", protein) #FUNC TEST

samples = ["S0"+str(i) for i in range(1,10)] + ["S10"] 
species = ["ARATH", "CAEEL", "HUMAN"]

df = df_triq #variable
df = df_spec #variable

sample = samples[0] #Variable
specie = "HUMAN" #variable
protein = "Q8TBA6_HUMAN" #VAR

protein = proteins[0] #iteration variable (if the protein is diff exp)
    
get_protein_abundance(triq, "S01", "HUMAN", protein) #FUNC TEST


#for protein in proteins:
#    for sample1 in samples:
#        abundance1 = get_protein_abundance(df, sample1, specie, protein)
#        for sample2 in samples:
#            abundance2 = get_protein_abundance(df, sample2, specie, protein)
            
#log2FC = np.log2(abundance1.values) - np.log2(abundance2.values)            
import time
start = time.time()

sample1 = "S02" #Var
sample2 = "S06" #Var
specie = "HUMAN" #Var
protein = proteins[0] #Var

log2FC_array = []
for protein in proteins:
    abundance_sample1 = get_protein_abundance(df, sample1, specie, protein)
    abundance_sample2 = get_protein_abundance(df, sample2, specie, protein)
    log2FC = np.log2(abundance_sample2) - np.log2(abundance_sample1)
    log2FC_array.append(log2FC)

end = time.time()
print(end-start)

# Process takes 15min


# Code to see that same proteins can have different indices
df_specie = df[df["specie"] == "HUMAN"]
df_specie_sample = df_specie[df_specie["sample"] == "S06"]
df_specie_sample = df_specie[df_specie["sample"] == "S02"]
df_specie_sample_protein = df_specie_sample[df_specie_sample["protein"] == protein]

from math import sqrt
test= []
for i in range(10):
    test.append(sqrt(i**2))

from math import sqrt
from joblib import Parallel, delayed
test2 = Parallel(n_jobs=2)(delayed(sqrt)(i ** 2) for i in range(10))

g = "g"
def get_val(x,g):
    return str(x)+str("_tetst")+str(g)


# write about this timing thingy here... wierd...
start=time.time()
test3 = Parallel(n_jobs=1)(delayed(get_val)(i,g) for i in ["a","b","c","d","e","b","c","d","e","b","c","d","e"])
end=time.time()
print(end-start)


# delayed(func)(parameter to func)

#### Make a function for abundance computations

def parallel_get_protein(df, protein, sample1, sample2, specie):
    abundance_sample1 = get_protein_abundance(df, sample1, specie, protein)
    abundance_sample2 = get_protein_abundance(df, sample2, specie, protein)
    log2FC = np.log2(abundance_sample2.values) - np.log2(abundance_sample1.values)
    return log2FC

import multiprocessing
num_cores = multiprocessing.cpu_count()


sample1 = "S02" #VAR
sample2 = "S06" #VAR
specie = "HUMAN" #VAR
df = spec #VAR




from log2fc_from_melt import get_log2FC_matrix 

df_log2FC = pd.DataFrame(log2FC_array, index=proteins, columns=runs).to_csv("spec_log2FC_S02_S06_HUMAN.csv", sep = "\t", index=False)
# Run this when we are at gym - Non-normalized but log2FC...
    
spec_log2FC_S02_S06_HUMAN.csv

print(1)
df_log2FC = get_log2FC_matrix(triq, "S02", "S06", "HUMAN")
df_log2FC.to_csv("triq_log2FC_S02_S06_HUMAN.csv", sep = "\t", index=False)
print(2)
df_log2FC = get_log2FC_matrix(triq, "S02", "S06", "CAEEL")
df_log2FC.to_csv("troq_log2FC_S02_S06_CAEEL.csv", sep = "\t", index=False)
print(3)
df_log2FC = get_log2FC_matrix(triq, "S02", "S06", "ARATH")
df_log2FC.to_csv("troq_log2FC_S02_S06_ARATH.csv", sep = "\t", index=False)
print(4)
df_log2FC = get_log2FC_matrix(spec, "S02", "S06", "CAEEL")
df_log2FC.to_csv("spec_log2FC_S02_S06_CAEEL.csv", sep = "\t", index=False)
print(5)
df_log2FC = get_log2FC_matrix(spec, "S02", "S06", "ARATH")
df_log2FC.to_csv("spec_log2FC_S02_S06_ARATH.csv", sep = "\t", index=False)

print(6)
df_log2FC = get_log2FC_matrix(triq, "S03", "S04", "HUMAN")
df_log2FC.to_csv("triq_log2FC_S03_S04_HUMAN.csv", sep = "\t", index=False)
print(7)
df_log2FC = get_log2FC_matrix(triq, "S03", "S04", "CAEEL")
df_log2FC.to_csv("triq_log2FC_S03_S04_CAEEL.csv", sep = "\t", index=False)
print(8)
df_log2FC = get_log2FC_matrix(triq, "S03", "S04", "ARATH")
df_log2FC.to_csv("triq_log2FC_S03_S04_ARATH.csv", sep = "\t", index=False)
print(9)
df_log2FC = get_log2FC_matrix(spec, "S03", "S04", "HUMAN")
df_log2FC.to_csv("spec_log2FC_S03_S04_HUMAN.csv", sep = "\t", index=False)
print(10)
df_log2FC = get_log2FC_matrix(spec, "S03", "S04", "CAEEL")
df_log2FC.to_csv("spec_log2FC_S03_S04_CAEEL.csv", sep = "\t", index=False)
print(11)
df_log2FC = get_log2FC_matrix(spec, "S03", "S04", "ARATH")
df_log2FC.to_csv("spec_log2FC_S03_S04_ARATH.csv", sep = "\t", index=False)

print(12)
df_log2FC = get_log2FC_matrix(triq, "S05", "S08", "HUMAN")
df_log2FC.to_csv("triq_log2FC_S05_S08_HUMAN.csv", sep = "\t", index=False)
print(13)
df_log2FC = get_log2FC_matrix(triq, "S05", "S08", "CAEEL")
df_log2FC.to_csv("triq_log2FC_S05_S08_CAEEL.csv", sep = "\t", index=False)
print(14)
df_log2FC = get_log2FC_matrix(triq, "S05", "S08", "ARATH")
df_log2FC.to_csv("triq_log2FC_S05_S08_ARATH.csv", sep = "\t", index=False)
print(15)
df_log2FC = get_log2FC_matrix(spec, "S05", "S08", "HUMAN")
df_log2FC.to_csv("spec_log2FC_S05_S08_HUMAN.csv", sep = "\t", index=False)
print(16)
df_log2FC = get_log2FC_matrix(spec, "S05", "S08", "CAEEL")
df_log2FC.to_csv("spec_log2FC_S05_S08_CAEEL.csv", sep = "\t", index=False)
print(17)
df_log2FC = get_log2FC_matrix(spec, "S05", "S08", "ARATH")
df_log2FC.to_csv("spec_log2FC_S05_S08_ARATH.csv", sep = "\t", index=False)
print(18)
df_log2FC = get_log2FC_matrix(triq, "S04", "S09", "HUMAN")
df_log2FC.to_csv("triq_log2FC_S04_S09_HUMAN.csv", sep = "\t", index=False)
print(19)
df_log2FC = get_log2FC_matrix(triq, "S04", "S09", "CAEEL")
df_log2FC.to_csv("triq_log2FC_S04_S09_CAEEL.csv", sep = "\t", index=False)
print(20)
df_log2FC = get_log2FC_matrix(triq, "S04", "S09", "ARATH")
df_log2FC.to_csv("triq_log2FC_S04_S09_ARATH.csv", sep = "\t", index=False)
print(21)
df_log2FC = get_log2FC_matrix(spec, "S04", "S09", "HUMAN")
df_log2FC.to_csv("spec_log2FC_S04_S09_HUMAN.csv", sep = "\t", index=False)
print(22)
df_log2FC = get_log2FC_matrix(spec, "S04", "S09", "CAEEL")
df_log2FC.to_csv("spec_log2FC_S04_S09_CAEEL.csv", sep = "\t", index=False)
print(23)
df_log2FC = get_log2FC_matrix(spec, "S04", "S09", "ARATH")
df_log2FC.to_csv("spec_log2FC_S04_S09_ARATH.csv", sep = "\t", index=False)
print("DONE!)

#try the pool approach

start=time.time()
log2FC_array_ = Parallel(n_jobs=num_cores)(delayed(parallel_get_protein)(df, protein, sample1, sample2, specie) for protein in proteins)
end=time.time()
print(time.time())






# compute differential expression of the log2FC dataframes.

ARATH_FC_matrix, CAEEL_FC_matrix, HUMAN_FC_matrix = get_log2FC_ratio_matrices()
ARATH_FC_matrix
CAEEL_FC_matrix



HUMAN_FC_matrix["S06"]["S02"]
spec_log2FC = pd.read_csv("spec_log2FC_S02_S06_HUMAN.csv", sep = "\t")
triq_log2FC = pd.read_csv("triq_log2FC_S02_S06_HUMAN.csv", sep = "\t")

method = "triq" # VAR
sample1 = "S02" # VAR
sample2 = "S06" # VAR
specie = "HUMAN" # VAR
ratio = 0.8 # VAR - ratio for what is considered diff. exp.
two_sided = False # VAR - is this needed? I will skip this for now



method = "triq" # VAR
sample1 = "S02" # VAR
sample2 = "S06" # VAR
specie = "HUMAN" # VAR
ratio = 0.8 # VAR - ratio for what is considered diff. exp.
two_sided = False # VAR - is this needed? I will skip this for now

from log2fc_from_melt_computation import get_differentially_expressed_proteins_from_log2FC_df


np.sum(get_differentially_expressed_proteins_from_log2FC_df("triq", "S02", "S06", "HUMAN", 0.8))
np.sum(get_differentially_expressed_proteins_from_log2FC_df("spec", "S02", "S06", "HUMAN", 0.8))

np.sum(get_differentially_expressed_proteins_from_log2FC_df("triq", "S02", "S06", "ARATH", 0.8))
np.sum(get_differentially_expressed_proteins_from_log2FC_df("spec", "S02", "S06", "ARATH", 0.8))

np.sum(get_differentially_expressed_proteins_from_log2FC_df("triq", "S02", "S06", "CAEEL", 0.8))
np.sum(get_differentially_expressed_proteins_from_log2FC_df("spec", "S02", "S06", "CAEEL", 0.8))

# Compute similar thing for proteinXvsY.csv - the posterior diff exp.






































