#!/usr/bin/env python3

import os.path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dataDir = "../intermediate"
resDir = "../result"
methodsDir = ["DID", "PS"]
eCol = "Actual Error Rate"
negCol = "Negative"

def countProteins(df, pCol, methodName):
    negative = df[pCol].str.contains("_HUMAN").astype(int)
    df[negCol] = negative.cumsum()
    df["pos"] = (1-negative).cumsum()
    df[eCol] = df[negCol]/(df["pos"]+df[negCol])
    df["method"] = methodName
    return df

def readMethods(resPath):
    msstats_raw = pd.read_csv(os.path.join(resPath, "msstat_output.csv"),
                                usecols=["Protein","pvalue","adj.pvalue"]) 
    msstats_raw.sort_values("pvalue", inplace=True)
    msstats_raw.rename(columns={'adj.pvalue' : 'fdr'}, inplace=True)
    msq_raw = pd.read_csv(os.path.join(resPath, "msqrob2_results.tsv"),
                                names=["Protein","logFC","se","df","t","pval","adjPval"], skiprows=1, usecols=["Protein","pval","adjPval"]) 
    msq_raw.dropna(inplace=True)                            
    msq_raw.sort_values("pval", inplace=True)
    msq_raw.rename(columns={'adjPval' : 'fdr', "" : "Protein"}, inplace=True)
    names = ["q_value","posterior_error_prob","protein"] + [str(a) for a in range(80)]
    triq_raw = pd.read_csv(os.path.join(resPath, "triqler_results/fc_0.96"),
                                usecols=[0,2], sep='\t', skiprows=1, names=names) 
    triq_raw.drop(triq_raw[triq_raw["protein"].str.contains("DECOY_")].index, inplace=True)
    triq_raw.rename(columns={'q_value' : 'fdr'}, inplace=True)
    return triq_raw, msstats_raw, msq_raw

for mDir in methodsDir:
    print(f"Processing data from {mDir}")
    subd = os.path.join(dataDir, mDir)
    outd = os.path.join(resDir, mDir)
    raws = readMethods(subd)
    methodNames = ["Triqler", "MSStats", "MSqRobSum"]
    counted = [countProteins(raw, p_name, method) for raw, p_name, method in zip(raws,["protein","Protein","Protein"], methodNames)]
    pn = pd.concat([df[["pos",negCol,"method"]] for df in counted], ignore_index=True, sort=False)
    pn.drop(pn[pn[negCol]>100].index, inplace=True)
    er_fdr = pd.concat([df[[eCol, "fdr", "method"]] for df in counted], ignore_index=True, sort=False)
    er_fdr.drop(er_fdr[er_fdr[eCol]>0.2].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr["fdr"].diff(periods=-1)>0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr[eCol].diff()<0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    pn.drop_duplicates(subset=['pos', 'method'], keep='last', inplace=True)
    pn.drop_duplicates(subset=[negCol, 'method'], keep='last', inplace=True)
    #print(pn.to_string())
    #print(er_fdr.to_string())
    sns.lineplot(data=pn, x=negCol, y="pos", hue="method", ci=None)
    plt.xlim(0.,50.)
    plt.show()
    plt.savefig(os.path.join(outd, "pos_neg.png"))
    plt.clf()
    sns.lineplot(data=er_fdr, x=eCol, y="fdr", hue="method", ci=None)
    plt.xlim(0.,0.1)
    plt.ylim(0.,0.3)
    plt.savefig(os.path.join(outd, "er_fdr.png"))
    plt.show()
    plt.clf()

