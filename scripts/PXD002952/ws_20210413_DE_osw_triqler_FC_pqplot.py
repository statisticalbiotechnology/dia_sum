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

from scipy import stats

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
    
    
    
    
    
#################
# OSW DE ########
#################
    
# read in df as above
    

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


df = pd.read_csv(file, sep = "\t")
df = df[df.decoy != 1]
df = df[df.m_score < 0.01] # filter away crap, so all values should be good... we take average of top3 here
df["experiment_id"] = df["filename"].map(experiment_id_mapper)
df["sample_id"] = df["filename"].map(sample_id_mapper)
df = df[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]
df["specie"] = df.ProteinName.map(specie_mapper)


def read_in_and_filter(filename, m_score_treshold = 0.01):    
    df = pd.read_csv(filename, sep = "\t")
    df = df[df.decoy != 1]
    df = df[df.m_score < 0.01] # filter away crap, so all values should be good... we take average of top3 here
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    df = df[["experiment_id", "sample_id", "Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName"]]
    df["specie"] = df.ProteinName.map(specie_mapper)
    return df
    
# highest scoring vs top highest intenstiy peptide?
# We take top3 intensity because we filter away so much with m_score < 0.01, all values should be good.

import time

def get_protein_quantity_sample_df(df):
    protein_array = []
    quantity_array = []
    m_score_array = []
    start = time.time()
    print("Peptide selection by top m_score")
    for peptide in df.FullPeptideName.unique():
        idx_best_score_peptide = df[df.FullPeptideName == peptide].m_score.sort_values(ascending = True).index[0]
        protQuant = df[df.index == idx_best_score_peptide].Intensity.values[0]
        protein =  df[df.index == idx_best_score_peptide].ProteinName.values[0]
        m_score = df[df.index == idx_best_score_peptide].m_score.values[0]
        protein_array.append(protein)
        quantity_array.append(protQuant)
        m_score_array.append(m_score)
    end = time.time()
    print(end-start)
    
    experiment_id = df.experiment_id.unique()[0]
    sample_id = df.sample_id.unique()[0]
    df_proteinQuant_with_duplicate_proteins = pd.DataFrame(np.array([protein_array, quantity_array, m_score_array]).T,
                                                           columns = ["protein", "intensity", "m_score"])
    
    # select protein with highest m_Score
    protein_array = []
    quantity_array = []
    start = time.time()
    print("Protein selection by top m_score")
    for protein in df_proteinQuant_with_duplicate_proteins.protein.unique():
        idx_best_score = df_proteinQuant_with_duplicate_proteins[df_proteinQuant_with_duplicate_proteins.protein == protein].m_score.sort_values(ascending=True).index[0]
        protQuant = df_proteinQuant_with_duplicate_proteins[df_proteinQuant_with_duplicate_proteins.index == idx_best_score].intensity.values[0]
        protein_array.append(protein)
        quantity_array.append(protQuant)
    end = time.time()
    print(end-start)
    df_protQuant = pd.DataFrame(np.array([quantity_array]).T, index = protein_array, columns = [experiment_id + "_" + sample_id])
    return df_protQuant


        
df = read_in_and_filter("HYE124_TTOF6600_32fix_lgillet_I150211_002-Pedro_-_Sample_1_-_SW32_-_Repl1.mzML.tsv", m_score_treshold = 0.01)
df1 = get_protein_quantity_sample_df(df)

df = read_in_and_filter("HYE124_TTOF6600_32fix_lgillet_I150211_003-Pedro_-_Sample_2_-_SW32_-_Repl1.mzML.tsv", m_score_treshold = 0.01)
df2 = get_protein_quantity_sample_df(df)

pd.concat([df1, df2], axis = 1)

df_protQuants = pd.DataFrame()
for file in os.listdir():
    if file[-4:] == ".tsv":
        df = read_in_and_filter(file, m_score_treshold=0.01)
        df = get_protein_quantity_sample_df(df)
        df_protQuants = pd.concat([df_protQuants, df], axis = 1)
df_protQuants.to_csv("protQuant.csv", sep ="\t")

def col_to_mIdx(df):
    def map_values_to_value_list(value_list, values):
        """
        e.g. test below lists to see results.
        value_list = ["a", "b", "c"]
        values = ["c","c","a"]
        """
        return [value_list.index(x) for x in values]
    exp_mapper = lambda x : x.split("_")[0]
    sample_mapper = lambda x : x.split("_")[1]
    
    exp_col = df.columns.map(exp_mapper)
    exp_list = list(exp_col.unique())
    exp_code = map_values_to_value_list(exp_list, exp_col)
    
    sample_col = df.columns.map(sample_mapper)
    sample_list = list(sample_col.unique())
    sample_code = map_values_to_value_list(sample_list, sample_col)

    midx = pd.MultiIndex(levels=[sample_list, exp_list], 
                         codes=[sample_code, exp_code],
                         names=["sample", "experiment"])
    return midx

def midx_df(df):
    midx = col_to_mIdx(df)
    df.columns = midx
    return df

df_protQuants = midx_df(df_protQuants)

def get_sample_values(df_protQuants, sample):
    df_sample = df_protQuants.iloc[:,df_protQuants.columns.get_level_values("sample") == sample]
    df_sample = df_sample[df_sample.isna().sum(axis = 1) < 1+max_nan] #
    df_sample = df_sample.astype(float)
    return df_sample

df_sample1 = get_sample_values(df_protQuants, sample = "1")
df_sample2 = get_sample_values(df_protQuants, sample = "2")

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

intersecting_proteins = intersection(list(df_sample1.index), list(df_sample2.index)) #This is the same as droping nan and might be an unfair comparison.

df_sample1 = df_sample1[df_sample1.index.isin(intersecting_proteins)]
df_sample2 = df_sample2[df_sample2.index.isin(intersecting_proteins)]
tStat  = stats.ttest_ind(df_sample1, df_sample2, axis = 1)[0]
pVals = stats.ttest_ind(df_sample1, df_sample2, axis = 1)[1]
df_log2FC = pd.DataFrame(np.array([np.log2(df_sample2.mean(axis=1))-np.log2(df_sample1.mean(axis=1))]).T, columns = ["log2FC"])
df_log2FC["p-value"] = pVals
df_log2FC["protein"] = intersecting_proteins
specie_mapper = lambda x : x.split("_")[-1]
df_log2FC["specie"] = df_log2FC.protein.map(specie_mapper)
df_log2FC = df_log2FC.set_index("protein")


"""
de_array = []

for fc_tresh in np.arange(0,2,0.05):
    df_log2FC["p-value"]
    n_de_up = (df_log2FC[df_log2FC.specie == "HUMAN"].log2FC > fc_tresh).sum()
    n_de_down = (df_log2FC[df_log2FC.specie == "HUMAN"].log2FC < -fc_tresh).sum()
    de_array.append(n_de_up + n_de_down)

specie_list = ["HUMAN", "ECOLI", "YEAS8"]

de_human = []
de_yeast = []
de_ecoli = []
de_file = []

"""



###############
###############



fc_treshold = 1.0
pVal_treshold = 0.05

def get_log2FC_with_significance(df_log2FC, fc_treshold, pVal_treshold):
    def significance(log2fc, pVal, fc_treshold = 1.0, pVal_treshold = 0.05):
        if log2fc > fc_treshold:
            if -np.log10(pVal) > -np.log10(pVal_treshold):
                return "up-regulated"
            else:
                return "non-significant"
        if log2fc < -fc_treshold:
            if -np.log10(pVal) > -np.log10(pVal_treshold):
                return "down-regulated"
            else:
                return "non-significant"
        else:
            return "non-significant"   
    df_log2FC["sign"] = df_log2FC.apply(lambda x: significance(x["log2FC"], x["p-value"], fc_treshold = fc_treshold, pVal_treshold = pVal_treshold), axis = 1)
    return df_log2FC

specie = "HUMAN"
fc = 1.0
pv = 0.05


def get_DE(df_log2FC, specie, fc, pv):
    df_sign = get_log2FC_with_significance(df_log2FC[df_log2FC.specie == specie], fc_treshold = fc, pVal_treshold = pv)
    de = (df_sign.sign != "non-significant").sum()
    return de

pv = 0.05
specie = "HUMAN"
def get_de_series(df_log2FC, specie, pv = 0.05):
    de_array = []
    for fc in np.arange(0,2.05,0.05):
        de = get_DE(df_log2FC, specie, fc, pv)
        de_array.append(de)
    return de_array


de_human = get_de_series(df_log2FC, "HUMAN", pv = 0.05)
de_ecoli = get_de_series(df_log2FC, "ECOLI", pv = 0.05)
de_yeast = get_de_series(df_log2FC, "YEAS8", pv = 0.05)

idx = pd.Series(np.arange(0,2.05,0.05))
idx.name = "fc_eval"
osw_DE = pd.DataFrame(np.array([de_human, de_ecoli, de_yeast]).T, columns = ["HUMAN", "ECOLI", "YEAST"], index = idx)

import seaborn as sns
import matplotlib.pyplot as plt
osw_DE = osw_DE.astype(float)
osw_DE.plot(grid=True)
plt.title("OSW (TOP3) n-differential expression")
###################
###################
## TRIQLER DE #####
###################
###################
    
    
import os 

import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

from triqler_output_to_df import  parse_triqler
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_results")


def get_triqler_DE(fc_tresh = 0.05): #WHGOOPS NOTE fc_tresh here is actually q-value treshold
    """
    Should be run in triqler results folder with files named fc_<fc>
    """
    specie_mapper = lambda x : x.split("_")[-1]
    specie_list = ["HUMAN", "ECOLI", "YEAS8"]
    
    de_human = []
    de_yeast = []
    de_ecoli = []
    de_file = []
    
    for file in sorted(os.listdir()):
        df = parse_triqler(file)
        df["specie"] = df.protein.map(specie_mapper)
        for specie in specie_list:
            df_specie = df[df.specie == specie]
            de_protein = len(df_specie[df_specie.q_value < fc_tresh])
            if specie == "HUMAN":
                de_human.append(de_protein)
            elif specie == "ECOLI":
                de_ecoli.append(de_protein)
            elif specie == "YEAS8":
                de_yeast.append(de_protein)
            else:
                print("ERROR")
        de_file.append(file)
    
    
    fc_mapper = lambda x : float(x.split("_")[-1])
    df_DE = pd.DataFrame(np.array([de_file, de_human, de_yeast, de_ecoli]), index = ["file", "HUMAN", "YEAST", "ECOLI"]).T
    df_DE["fc_eval"]= df_DE.file.map(fc_mapper)
    return df_DE

df_triq_DE = get_triqler_DE(fc_tresh = 0.05)


import seaborn as sns
import matplotlib.pyplot as plt
df_triq_DE = df_triq_DE.drop("file", axis = 1)
df_triq_DE = df_triq_DE.set_index("fc_eval")
df_triq_DE = df_triq_DE.astype(float)
df_triq_DE.plot(grid=True)
plt.title("Triqler n-differential expression")





