#!/usr/bin/env python3

import os.path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_context("talk")

dataDir = "../intermediate"
resDir = "../result"
methodsDir = ["DID", "PS"]
eCol = "Actual Error Rate"
negCol = "Differential HeLa"
posCol = "Differential non-HeLa"

def countProteins(df, pCol, methodName):
    negative = df[pCol].str.contains("_HUMAN").astype(int)
    df[negCol] = negative.cumsum()
    df[posCol] = (1-negative).cumsum()
    df[eCol] = df[negCol]/(df[posCol]+df[negCol])
    df["method"] = methodName
    return df

def readMethods(resPath):
    msstats_raw = pd.read_csv(os.path.join(resPath, "msstat_output.csv"),
                                usecols=["Protein","pvalue","adj.pvalue"]) 
    msstats_raw.sort_values("pvalue", inplace=True)
    msstats_raw.rename(columns={'adj.pvalue' : "FDR"}, inplace=True)
    msq_raw = pd.read_csv(os.path.join(resPath, "msqrob2_results.tsv"),
                                names=["Protein","logFC","se","df","t","pval","adjPval"], skiprows=1, usecols=["Protein","pval","adjPval"]) 
    msq_raw.dropna(inplace=True)                            
    msq_raw.sort_values("pval", inplace=True)
    msq_raw.rename(columns={'adjPval' : "FDR", "" : "Protein"}, inplace=True)
    names = ["q_value","posterior_error_prob","protein"] + [str(a) for a in range(80)]
    triq_raw = pd.read_csv(os.path.join(resPath, "triqler_results/fc_0.96"),
                                usecols=[0,2], sep='\t', skiprows=1, names=names) 
    triq_raw.drop(triq_raw[triq_raw["protein"].str.contains("DECOY_")].index, inplace=True)
    triq_raw.rename(columns={'q_value' : "FDR"}, inplace=True)
    return triq_raw, msstats_raw, msq_raw

for mDir in methodsDir:
    print(f"Processing data from {mDir}")
    subd = os.path.join(dataDir, mDir)
    outd = os.path.join(resDir, mDir)
    raws = readMethods(subd)
    methodNames = ["Triqler", "MSStats", "MSqRobSum"]
    counted = [countProteins(raw, p_name, method) for raw, p_name, method in zip(raws,["protein","Protein","Protein"], methodNames)]
    pn = pd.concat([df[[posCol,negCol,"method"]] for df in counted], ignore_index=True, sort=False)
    pn.drop(pn[pn[negCol]>100].index, inplace=True)
    er_fdr = pd.concat([df[[eCol, "FDR", "method"]] for df in counted], ignore_index=True, sort=False)
    er_fdr.drop(er_fdr[er_fdr[eCol]>0.2].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr["FDR"].diff(periods=-1)>0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr[eCol].diff()<0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    pn.drop_duplicates(subset=[posCol, 'method'], keep='last', inplace=True)
    pn.drop_duplicates(subset=[negCol, 'method'], keep='last', inplace=True)
    #print(pn.to_string())
    #print(er_fdr.to_string())
    g = sns.lineplot(data=pn, x=negCol, y=posCol, hue="method", ci=None)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles=handles[1:], labels=labels[1:])
    g.get_legend().set_title(None)
    plt.xlim(0.,50.)
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(outd, "pos_neg.png"))
    g = sns.lineplot(data=er_fdr, x=eCol, y="FDR", hue="method", ci=None)
    g.get_legend().set_title(None)
    plt.xlim(0.,0.1); plt.ylim(0.,0.3)
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(outd, "er_fdr.png"))

