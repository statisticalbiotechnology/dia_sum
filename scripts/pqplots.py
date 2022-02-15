#!/usr/bin/env python3

import os.path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dataDir = "../intermediate"
resDir = "../result"
methodsDir = ["DID", "PS"]

def countProteins(df, pCol, methodName):
    print(df)
    negative = df[pCol].str.contains("_HUMAN").astype(int)
    df["neg"] = negative.cumsum()
    df["pos"] = (1-negative).cumsum()
    df["error_rate"] = df["neg"]/df["pos"]
    df["method"] = methodName
    return df

def readMethods(resPath):
    msstats_raw = pd.read_csv(os.path.join(resPath, "msstat_output.csv"),
                                usecols=["Protein","pvalue","adj.pvalue"]) 
    msstats_raw.sort_values("pvalue", inplace=True)
    msq_raw = pd.read_csv(os.path.join(resPath, "msstat_output.csv"),
                                usecols=["Protein","pvalue","adj.pvalue"]) 
    msq_raw.dropna(inplace=True)                            
    msq_raw.sort_values("pvalue", inplace=True)
    names = ["q_value","posterior_error_prob","protein"] + [str(a) for a in range(80)]
    triq_raw = pd.read_csv(os.path.join(resPath, "triqler_results/fc_0.96"),
                                usecols=[0,2], sep='\t', names=names) 
    # triq_raw = pd.read_csv(os.path.join(resPath, "triqler_results/fc_0.96"),
    #                             usecols=["protein","q_value"], sep='\t') 
    return triq_raw, msstats_raw, msq_raw

for mDir in methodsDir:
    subd = os.path.join(dataDir, mDir)
    outd = os.path.join(resDir, mDir)
    raws = readMethods(subd)
    methodNames = ["triq","mss","msq"]
    counted = [countProteins(raw, p_name, method) for raw, p_name, method in zip(raws,["protein","Protein","Protein"], methodNames)]
    pn = pd.concat([df[["pos","neg","method"]] for df in counted], ignore_index=True, sort=False)
    er_fdr = pd.concat([df[["error_rate", qName, "method"]] for df, qName in zip(counted, ["q_value", "adj.pvalue", "adj.pvalue"])], ignore_index=True, sort=False)
    er_fdr["fdr"] = er_fdr[["adj.pvalue","q_value"]].max() # One column is NA, select the other
    sns.lineplot(data=pn, x="neg", y="pos", hue="method")
    plt.savefig(os.path.join(outd, "pos_neg.png"))
    sns.lineplot(data=er_fdr, x="error_rate", y="fdr", hue="method")
    plt.savefig(os.path.join(outd, "er_fdr.png"))

