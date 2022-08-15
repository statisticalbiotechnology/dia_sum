#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 17:51:36 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import os 

os.chdir("/hdd_14T/data/pxd021197/ftp.ebi.ac.uk/pride-archive/2021/07/PXD021197/DIA")

dfs = []
for i in os.listdir():
    if i[-4:] == ".tsv":
        if i[:2] == "CK":
            df = pd.read_csv(i, sep = "\t")
            dfs.append(df)

df = pd.concat(dfs)
df.to_csv("result.tsv", sep = "\t", index = False)



def control_condition_rename(x):
    if x in ["Healthy control 1", "Healthy control 2", "Healthy control 3", "Healthy control 4", "Bacterial pneumonia control"]:
        return "control"
    else:
        return x

def control_case_rename(x):
    if x not in ["TP1", "control"]:
        return "Case"
    else:
        return "Control"
    
def convert_pxd021197_to_triqler(filename = "plasma.aligned.tsv", mapper_file = "PXD021197_DIA_DDA.xlsx", output = "output.tsv"):
    mapper = pd.read_excel(mapper_file)
    df = pd.read_csv(filename, sep = "\t", usecols = ["Charge", "m_score", "Intensity", "FullPeptideName", "ProteinName", "filename"])
    df["mapper"] = df.filename.map(lambda x:x.split(".")[0])
    mapper.SampleName
    len(mapper.Patient.unique())
    mapper.columns
    patient_map=mapper[["FileName", "Patient"]].set_index("FileName").to_dict()["Patient"]
    run_map=mapper[["FileName", "SampleName"]].set_index("FileName").to_dict()["SampleName"]
    df["patient"] = df.mapper.map(patient_map)
    #df.condition = df.condition.map(control_condition_rename)
    # case control rename
#    df.condition = df.condition.map(control_case_rename)
    df["run"] = df.mapper.map(run_map)
    df["condition"] = df.run.map(lambda x:x.split("_")[2].split("+")[0].split("-")[0])
    #df["condition"] = df.run.map(lambda x:x.split("_")[2]) # every day split
    df.rename({"Charge":"charge", "m_score":"searchScore", "Intensity":"intensity", "FullPeptideName":"peptide", "ProteinName":"proteins"},axis=1, inplace = True)
    df.columns
    df = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    df["searchScore"] = -np.log10(df["searchScore"])

    df.intensity = np.log(df.intensity).copy()
    return df

def remove_sputum(res):
    res = res[~res.run.str.contains("SpPlug")].copy()
    return res

def remove_healhty(res):
    res = res[~res.run.str.contains("Hea")].copy()
    return res

def remove_patpne(res):
    res = res[-res.run.str.contains("PatPne_Plasma_NA")].copy()
    return res

def remove_single_run_per_condition_day_split(res):
    res = res[res.run != "TP3_Plasma_Post+8"].copy()
    res = res[res.run != "TP5_Plasma_Pre-4"].copy()
    return res





df = df[~df.run.str.contains("TP3_Plasma_Post+8")].copy()

df[df.run.str.contains("TP3_Plasma_Post+8")]
df[df.condition.str.contains´("Post+8")]
df.condition.unique()
df.run.unique()


df = remove_sputum(df)
df = remove_healhty(df)
df = remove_patpne(df)

# remove two runs
df = remove_single_run_per_condition_day_split(df)
df.groupby("condition").run.unique()


res.condition.unique()

os.chdir("/home/ptruong/Downloads/pxd021197")

filename = "plasma.aligned.tsv"
filename = "result.tsv"
mapper_file = "PXD021197_DIA_DDA.xlsx"

df.decoy.sum()

res_runs = df.run.unique()
df[df.run.isin(res_runs)].proteins.str.contains("DECOY").sum()

df
df[~df.run.str.contains("TP3_Plasma_Post+8")].run.unique()
df.run.unique()
df[df.run.isin(res_runs)].to_csv("triqler_input_pre_post_reverse_mscore_log.tsv", sep = "\t", index = False)

df_t = pd.read_csv("triqler_input.tsv", sep = "\t")
df_t.intensity

df.condition = df.condition.map(lambda x:x.split("+")[0])


df_pre_post = df[df.condition.isin(["Post+7" ,"Post+6" ,"Post+5", "Post+4", "Post+3", "Pre-1", "Pre-2", "Pre-3"])]
df_pre_post.condition = df_pre_post.condition.map(lambda x:x.split("+")[0])
df_pre_post.to_csv("triqler_input_pre_post_reverse_mscore_trunc.tsv", sep = "\t", index = False)

