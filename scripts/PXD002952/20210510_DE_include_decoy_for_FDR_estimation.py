#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:13:46 2021

@author: ptruong
"""



import pandas as pd
import numpy as np

from scipy import stats

import os 

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")

from q_value import qvalues
from triqler_output_to_df import  parse_triqler
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/full_ts_v")


# filename has different formatting, we need to change number or implement regex.
experiment_id_mapper = lambda x: x.split("_")[5]
sample_id_mapper = lambda x: x.split("_")[8] #hye124 
specie_mapper = lambda x: x.split("_")[-1]

def decoy_mapper(x):
    if x.split("_")[0] == "DECOY":
        return True
    else:
        return False
    
def read_in_and_filter(filename, m_score_treshold = 0.01):  
    print(filename)
    df = pd.read_csv(filename, sep = "\t")
    #df = df[df.decoy != 1]
    df = df[df.m_score < m_score_treshold] # filter away crap, so all values should be good... we take average of top3 here
    print(str(len(df)) + " significantly identified peptides at " + str(m_score_treshold) + " FDR-treshold.")
    print("")
    df["experiment_id"] = df["filename"].map(experiment_id_mapper)
    df["sample_id"] = df["filename"].map(sample_id_mapper)
    sample_id = df.sample_id[0]
    experiment_id = df.experiment_id[0]     
    def top3(df):
        df = (df.groupby('ProteinName')['Intensity'].apply(lambda x: x.nlargest(3).mean() if len(x.nlargest(3)) >= 1 else np.nan)
                  .reset_index())
        #print(df.isna().sum())
        return df
    df_reduced = df[["ProteinName", "Intensity"]]
    df_protein = top3(df_reduced)
    df = df_protein
    df["specie"] = df.ProteinName.map(specie_mapper)
    midx = pd.MultiIndex(levels = [[sample_id],[experiment_id]], codes = [[0],[0]], names = ["sample_id", "experiment_id"])
    df = df.set_index(["specie", "ProteinName"])
    df = pd.DataFrame(df.values, columns = midx, index = df.index)
    
    return df



# highest scoring vs top highest intenstiy peptide?
# We take top3 intensity because we filter away so much with m_score < 0.01, all values should be good.

dfs = []
for file in os.listdir():
    if file[-8:] == "mzML.tsv":
        dfs.append(read_in_and_filter(file, m_score_treshold=0.01))
        #print(len(df_part))
        #df = pd.concat([df, df_part],axis = 1)        
df = pd.concat(dfs, axis = 1)

tmp = dfs[0]
tmp
tmp["decoy"] = tmp.index.get_level_values("ProteinName").map(decoy_mapper)
tmp[tmp.decoy == True]["2"]["005-Pedro"].unique()
len(tmp[tmp.decoy == True]["2"]["005-Pedro"].unique())


A = df[df.iloc[:, df.columns.get_level_values("sample_id") == "1"].isna().sum(axis=1)<2]
A = A.iloc[:,A.columns.get_level_values("sample_id") == "1"]
#A = np.log2(A)
B = df[df.iloc[:, df.columns.get_level_values("sample_id") == "2"].isna().sum(axis=1)<2]
B = B.iloc[:,B.columns.get_level_values("sample_id") == "2"]
#B = np.log2(B)



#A["decoy"] = A.index.get_level_values("ProteinName").map(decoy_mapper)
#B["decoy"] = B.index.get_level_values("ProteinName").map(decoy_mapper)
#A[A.decoy == True]
#B[B.decoy == True]

# Find overlapping proteins
overlapping_proteins = list(set(A.index) & set(B.index))
A = A[A.index.isin(overlapping_proteins)]
B = B[B.index.isin(overlapping_proteins)]


import scipy.stats as stats

p_vals = stats.ttest_ind(A, B, axis = 1)[1]
p_vals = pd.DataFrame(p_vals, columns = ["p"])
p_vals["q"] = qvalues(p_vals)
p_vals = pd.DataFrame(p_vals.values, index = A.index, columns = ["p", "q"])
p_vals.sort_values("q",inplace =True)
p_vals = p_vals.astype(float)

A = np.log2(A.sum(axis=1))
B = np.log2(B.sum(axis=1))

A.name = "1"
B.name = "2"

df_final = pd.concat([A, B, p_vals], axis = 1)
df_final["log2(A,B)"] = df_final["1"] - df_final["2"]
df_final["decoy"] = df_final.index.get_level_values("ProteinName").map(decoy_mapper)

def pq_data(df_final, fc_treshold):
    # two-side treshold because triqler uses two side treshold.
    df_final_fc_gt = (df_final[df_final["log2(A,B)"] > fc_treshold])
    df_final_fc_lt = (df_final[df_final["log2(A,B)"] < -fc_treshold])
    df_final_fc = pd.concat([df_final_fc_gt, df_final_fc_lt])
    n_array = []
    for q in np.arange(0,0.101, 0.001):
        n = (df_final_fc["q"] <= q).sum()
        n_array.append(n)
    res = pd.DataFrame(n_array, index = np.arange(0,0.101, 0.001), columns = ["DE"])
    return res

# test
#pq_data(df_final[df_final.index.get_level_values("specie") == "HUMAN"], fc_treshold = 0.2)
#pq_data(df_final[df_final.index.get_level_values("specie") == "YEAS8"], fc_treshold = 0.2)
#pq_data(df_final[df_final.index.get_level_values("specie") == "ECOLI"], fc_treshold = 0.2)


import matplotlib.pyplot as plt

specie = "HUMAN"

def plot_pq(df_final, specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    for i in range(10):
        fc =i*0.2        
        df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)

        if i == 0:
            y_lim = df_top3_decoy.max()[0]
        #axs[row, col].plot(df_top3.index, df_top3.DE)
        axs[row, col].plot(df_top3_decoy.index, df_top3_decoy.DE)

        axs[row, col].set_ylim([0, y_lim])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("n - Differentially expressed.")
        col += 1
        if col == 5:
            row +=1
            col = 0
    plt.suptitle("OSW top3 " + specie)
    plt.show()
    
df_final[df_final.decoy == True]

plot_pq(df_final, specie = "HUMAN")
plot_pq(df_final, specie = "ECOLI")
plot_pq(df_final, specie = "YEAS8")


import time
from triqler_output_to_df import  parse_triqler

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_results_20210510")

def get_pq_data_triqler(q_val):
    fc_tresh = []
    n_hs = []
    n_ye = []
    n_ec = []
    
    n_hs_decoy = []
    n_ye_decoy = []
    n_ec_decoy = []
    
    for file in sorted(os.listdir()): 
        if file[:3] == "fc_":
            fc = float(file.split("_")[1])
            df_triq = parse_triqler(file)
            df_triq = df_triq[df_triq.q_value < q_val]
            df_triq["specie"] = df_triq.protein.map(specie_mapper)
            df_triq["decoy"] = df_triq.protein.map(decoy_mapper)
            n_hs.append((df_triq[df_triq["decoy"] == False]["specie"] == "HUMAN").sum())
            n_ye.append((df_triq[df_triq["decoy"] == False]["specie"] == "YEAS8").sum())
            n_ec.append((df_triq[df_triq["decoy"] == False]["specie"] == "ECOLI").sum())
            
            n_hs_decoy.append((df_triq[df_triq["decoy"] == True]["specie"] == "HUMAN").sum())
            n_ye_decoy.append((df_triq[df_triq["decoy"] == True]["specie"] == "YEAS8").sum())
            n_ec_decoy.append((df_triq[df_triq["decoy"] == True]["specie"] == "ECOLI").sum())
            
            fc_tresh.append(fc)
    
    df =pd.DataFrame(np.array([n_hs, n_ye, n_ec,
                               n_hs_decoy, n_ye_decoy, n_ec_decoy,]).T, 
    index = fc_tresh, columns = ["HUMAN", "YEAS8", "ECOLI",
                                 "HUMAN_decoy", "YEAS8_decoy", "ECOLI_decoy"])
    return df


#start = time.time()
def get_DE_for_fcs():
    qs = []
    dfs = []
    for q in np.arange(0,0.101, 0.001):
        qs.append(q)
        dfs.append(get_pq_data_triqler(q))
        #print(time.time()-start)
    end = time.time()
    #print(end-start)
    
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    res = []
    for fc in [round(i*0.2,2) for i in range(10)]:
        vals = []
        for df in dfs:
            val = list(df[df.index == fc].values[0])
            vals.append(val)    
        df_res = pd.DataFrame(vals, index = np.arange(0,0.101, 0.001), columns = ["HUMAN", "YEAS8", "ECOLI",
                               "HUMAN_decoy", "YEAS8_decoy", "ECOLI_decoy"])
        res.append(df_res)
    return res

res = get_DE_for_fcs()



#Rewrite this to function subplot
def plot_pq_specie(specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    #specie = "ECOLI"
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    for i in range(10):
        fc = fcs[i]
        res[i][specie].plot(ax = axs[row,col]) #triqler data plot
        df = pq_data(df_final[df_final.index.get_level_values("specie") == specie], fc_treshold = fc)
        #df_ms = msstats_pq_data(ms[ms.specie == specie], fc = fc)
        axs[row, col].plot(df.index, df.DE) # plot osw top3
        axs[row, col].plot(df_ms.index, df_ms.DE) # plot MSSTATS
        axs[row, col].legend(labels=["Triqler", "OSW Top3", "msStats"])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("n - Differentially expressed.")
        col+=1
        if col == 5:
            row+=1
            col=0
    plt.suptitle(f"Differentially expressed proteins {specie} ( True ratios 1:1 for human, 2:1 for yeast, and 1:4 for E.coli )")# + specie)
    plt.show()


def plot_pq_specie(specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    for i in range(10):
        fc =fcs[i]        
        df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        
        res[i][specie].plot(ax = axs[row,col], style = "b") #triqler data plot
        res[i][specie+"_decoy"].plot(ax = axs[row,col], style = "b--") #triqler data plot
        #if i == 0:
        #    y_lim = df_top3_decoy.max()[0]
        axs[row, col].plot(df_top3.index, df_top3.DE, "g")
        axs[row, col].plot(df_top3_decoy.index, df_top3_decoy.DE, "g--")
        
        #axs[row, col].set_ylim([0, y_lim])
        axs[row, col].legend(labels=["Triqler", "Triqler_decoy", "Top3", "Top3 decoy"])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("n - Differentially expressed.")
        col += 1
        if col == 5:
            row +=1
            col = 0
    plt.suptitle("OSW top3 " + specie)
    plt.show()

plot_pq_specie(specie = "HUMAN")
plot_pq_specie(specie = "YEAS8")
plot_pq_specie(specie = "ECOLI")




### MAKE CALIBRATION PLOTS
fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
fcs.sort()

i=3
specie = "YEAS8"

fc = fcs[i]
df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)

fig, axs = plt.subplots(2, 5)
#fig, ax = plt.subplots()
axs.plot([0, 1], [0, 1], "--", transform=ax.transAxes)
(df_top3_decoy/ df_top3).plot(ax = ax)
(res[i][specie+"_decoy"]/res[i][specie]).plot(style = "--", ax= ax)


def plot_calibration_specie(specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    for i in range(10):
        fc =fcs[i]        
        df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        
        axs[row, col].plot([0, 1], [0, 1], "--", transform=axs[row,col].transAxes)
        (df_top3_decoy/ df_top3).plot(ax = axs[row, col])
        (res[i][specie+"_decoy"]/res[i][specie]).plot(style = "-", ax=axs[row,col] )

        #axs[row, col].set_ylim([0, y_lim])
        axs[row, col].legend(labels=["", "top3", "triqler"])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("DE FP/TP.")
        col += 1
        if col == 5:
            row +=1
            col = 0
    plt.suptitle("Calibration plot " + specie + "(True ratios: Human 1:1, Yeast 2:1, Ecoli 1:4)")
    plt.show()


plot_calibration_specie(specie = "HUMAN")
plot_calibration_specie(specie = "YEAS8")
plot_calibration_specie(specie = "ECOLI")

############
# MS STATS #
############


def msstats_pq_data(ms, fc):
    # ms = pd.read_csv("msstats.csv")
    ms_gt = ms[ms.log2FC > fc]
    ms_lt = ms[ms.log2FC < -fc]
    ms_fc = pd.concat([ms_gt, ms_lt])
    n_array = []
    for q in np.arange(0,0.101, 0.001):
        n = (ms_fc["adj.pvalue"] <= q).sum()
        n_array.append(n)
    res = pd.DataFrame(n_array, index = np.arange(0,0.101, 0.001), columns = ["DE"])
    return res

# MSSTATS
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")

ms = pd.read_csv("msstats_20210511_unfiltered.csv")
ms["specie"] = ms.ProteinName.map(specie_mapper)
ms["decoy"] = ms.ProteinName.map(decoy_mapper)
len(ms[ms.decoy == False])/len(ms[ms.decoy == True])

test = pd.read_csv("data.transition.csv", sep = ",")
test = pd.read_csv("data.annotated.csv", sep = ",")
test = pd.read_csv("data.filtered.csv", sep = ",")
test.ProteinName.map(decoy_mapper).sum()



def msstats_pq_data(ms, fc):
    # ms = pd.read_csv("msstats.csv")
    ms_gt = ms[ms.log2FC > fc]
    ms_lt = ms[ms.log2FC < -fc]
    ms_fc = pd.concat([ms_gt, ms_lt])
    n_array = []
    for q in np.arange(0,0.101, 0.001):
        n = (ms_fc["adj.pvalue"] <= q).sum()
        n_array.append(n)
    res = pd.DataFrame(n_array, index = np.arange(0,0.101, 0.001), columns = ["DE"])
    return res

# MSSTATS
os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")

ms = pd.read_csv("msstat_DE_res_20210512.csv")
ms["specie"] = ms.Protein.map(specie_mapper)
ms["decoy"] = ms.Protein.map(decoy_mapper)


def plot_pq_specie(specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    for i in range(10):
        fc =fcs[i]        
        df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        
        df_ms = msstats_pq_data(ms[(ms.specie == specie) & (ms.decoy == False)], fc = fc)
        df_ms_decoy = msstats_pq_data(ms[(ms.specie == specie) & (ms.decoy == True)], fc = fc)

        res[i][specie].plot(ax = axs[row,col], style = "b") #triqler data plot
        res[i][specie+"_decoy"].plot(ax = axs[row,col], style = "b--") #triqler data plot
        #if i == 0:
        #    y_lim = df_top3_decoy.max()[0]
        axs[row, col].plot(df_top3.index, df_top3.DE, "g")
        axs[row, col].plot(df_top3_decoy.index, df_top3_decoy.DE, "g--")
        
        # msstat
        axs[row, col].plot(df_ms.index, df_ms.DE, "r")
        axs[row, col].plot(df_ms_decoy.index, df_ms_decoy.DE, "r--")
        
        #axs[row, col].set_ylim([0, y_lim])
        axs[row, col].legend(labels=["Triqler", "Triqler_decoy", "Top3", "Top3 decoy", "msStat", "msStat decoy"])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("n - Differentially expressed.")
        col += 1
        if col == 5:
            row +=1
            col = 0
    plt.suptitle(f"Differentially expressed proteins {specie} ( True ratios 1:1 for human, 2:1 for yeast, and 1:4 for E.coli )")
    plt.show()

plot_pq_specie(specie = "HUMAN")
plot_pq_specie(specie = "YEAS8")
plot_pq_specie(specie = "ECOLI")

ms[(ms.specie == specie) & (ms.decoy == False)]

ms[(ms.specie == specie) & (ms.decoy == True)]


def plot_calibration_specie(specie):
    fig, axs = plt.subplots(2, 5)
    row = 0
    col = 0
    fcs = [round(i*0.2,2) for i in range(9)] + [0.68]
    fcs.sort()
    for i in range(10):
        fc =fcs[i]        
        df_top3 = pq_data(df_final[(df_final.decoy == False) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        df_top3_decoy = pq_data(df_final[(df_final.decoy == True) & (df_final.index.get_level_values("specie") == specie)], fc_treshold = fc)
        
        df_ms = msstats_pq_data(ms[(ms.specie == specie) & (ms.decoy == False)], fc = fc)
        df_ms_decoy = msstats_pq_data(ms[(ms.specie == specie) & (ms.decoy == True)], fc = fc)
        
        axs[row, col].plot([0, 1], [0, 1], "--", transform=axs[row,col].transAxes)
        (df_top3_decoy/ df_top3).plot(ax = axs[row, col])
        (res[i][specie+"_decoy"]/res[i][specie]).plot(style = "-", ax=axs[row,col] )
        (df_ms_decoy / df_ms).plot(ax = axs[row, col])

        #axs[row, col].set_ylim([0, y_lim])
        axs[row, col].legend(labels=["", "top3", "triqler", "msstat"])
        axs[row, col].set_title(f"fc = {fc}")
        axs[row, col].set_xlabel("q-value")
        axs[row, col].set_ylabel("DE FP/TP.")
        col += 1
        if col == 5:
            row +=1
            col = 0
    plt.suptitle("Calibration plot " + specie + "(True ratios: Human 1:1, Yeast 2:1, Ecoli 1:4)")
    plt.show()


plot_calibration_specie(specie = "HUMAN")
plot_calibration_specie(specie = "YEAS8")
plot_calibration_specie(specie = "ECOLI")
