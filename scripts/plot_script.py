#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:28:18 2022

@author: ptruong
"""

import re
import os
import os.path
import numpy as np
import numpy.random as npr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
sns.set_context("talk")

dataDir = "../intermediate"
dataDir = "/hdd_14T/data/PXD002952/20210614_dataset/result_files_20220214"

resDir = "../result"
resDir = "/hdd_14T/data/PXD002952/20210614_dataset/result_files_20220214/result"
methodsDir = ["DID", "PS"]
eCol = "Actual Error Rate"
negCol = "Differential HeLa"
posCol = "Differential non-HeLa"

#os.chdir(dataDir)

def bootstrap(invec):
    """ Function for generating bootstrap sampled versions of a vector """
    idx = npr.randint(0, len(invec), len(invec))
    return [invec[i] for i in idx]

def estimatePi0(p, numBoot=100, numLambda=100, maxLambda=0.95):
    """
    Function for estimaring pi_0, i.e. the prior null probability
    for p values.
    Args:
        p (list(float)): The list of p values for which pi_0 should be estimated
        numBoot (int): The number of bootstrap rounds that should be made.
        numLambda (int): The number of lambda tresholds that should be evaluated.
        maxLambda (float): The upper bond of the range of lambda treshold.
    Returns:
        pi_0, a float containing the pi_0 estimate.
    """
    p.sort()
    n=len(p)
    lambdas=np.linspace(maxLambda/numLambda,maxLambda,numLambda)
    Wls=np.array([n-np.argmax(p>=l) for l in lambdas])
    pi0s=np.array([Wls[i] / (n * (1 - lambdas[i])) for i in range(numLambda)])
    minPi0=np.min(pi0s)
    mse = np.zeros(numLambda)
    for boot in range(numBoot):
        pBoot = bootstrap(p)
        pBoot.sort()
        WlsBoot =np.array([n-np.argmax(pBoot>=l) for l in lambdas])
        pi0sBoot =np.array([WlsBoot[i] / (n *(1 - lambdas[i])) for i in range(numLambda)])
        mse = mse + np.square(pi0sBoot-minPi0)
    minIx = np.argmin(mse)
    return pi0s[minIx]

def qvalues(pvalues,p_col = "p", q_col ="q", pi0 = 1.0):
    """
    Function for estimaring q values.
    Args:
        pvalues (DataFrame): A DataFrame that contain at least one column with pvalues.
        p_col (str): The name of the column that contain pvalues.
        q_col (str): The name of the column that that shall contain the estimated
                    q-values. The column will be created if not already existing.
        pi0 (float): The prior probability of the null hypothesis. If set to None, this is estimated from data.
                    Defaults to 1.0
    Returns:
        The modified DataFrame.
    """
    m = pvalues.shape[0] # The number of p-values
    pvalues.sort_values(p_col,inplace=True) # sort the pvalues in acending order
    if pi0 is None:
        pi0 = estimatePi0(list(pvalues[p_col].values))

    # calculate a FDR(t) as in Storey & Tibshirani
    num_p = 0.0
    for ix in pvalues.index:
        num_p += 1.0
        fdr = pi0*pvalues.loc[ix,p_col]*m/num_p
        pvalues.loc[ix,q_col] = fdr

    # calculate a q(p) as the minimal FDR(t)
    old_q=1.0
    for ix in reversed(list(pvalues.index)):
        q = min(old_q,pvalues.loc[ix,q_col])
        old_q = q
        pvalues.loc[ix,q_col] = q
    return pvalues

def oswGenerator(path):
    for file in os.listdir(path):
        yield pd.read_csv(os.path.join(path,file), sep='\t', usecols =["ProteinName", "Intensity","Sequence"]), file, "ProteinName", "Intensity"

def diannGenerator(path):
    spec_raw = pd.read_csv(path, sep='\t', usecols =["Run", "Protein.Ids", "Precursor.Normalised","Modified.Sequence"])
    for run, rdata in spec_raw.groupby("Run"):
        yield rdata, run, "Protein.Ids", "Precursor.Normalised"

def top3ByGenerator(generator):
    pattern = re.compile(r'(Sample_)(\d+)(.+Repl)(\d+)')
    serieses = []
    for run, name, protCol, intensityCol in generator:
        mo = pattern.search(name)
        sample, repl = mo[2],mo[4]
        run_name = f"{sample}_{repl}"
        p2i = {}
        for protein, peptides in run.groupby(protCol):
            if peptides.shape[0]>0:
                peptides.sort_values(intensityCol, inplace=True, ascending=False)
                top3 = peptides[intensityCol][:min(3,peptides.shape[0])].mean()
                p2i[protein] = top3
        serieses.append(pd.Series(p2i, name=run_name))
    df = pd.concat(serieses, axis=1, join="inner")
    df.index.name = "Protein"
    df.reset_index(level=0, inplace=True)
    df["pval"] = stats.ttest_ind(df[[col for col in df.columns if '1_' in col]], 
                               df[[col for col in df.columns if '2_' in col]], 
                               axis = 1)[1]
    return qvalues(df,"pval","FDR", None)


def countProteins(df, pCol, methodName):
    df.drop(df[df[pCol].str.contains("DECOY_")].index, inplace=True)
    negative = df[pCol].str.contains("_HUMAN").astype(int)
    df[negCol] = negative.cumsum()
    df[posCol] = (1-negative).cumsum()
    df[eCol] = df[negCol]/(df[posCol]+df[negCol])
    df["method"] = methodName
    return df


def readMethods(resPath):
    # reads the different file formats and place them in separate data frames
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
    triq_raw.rename(columns={'q_value':"FDR", "protein":"Protein"}, inplace=True)
    if "DID" in resPath:
        generator = oswGenerator(os.path.join(resPath,"osw_results"))
    else:
        generator = specialGenerator(os.path.join(resPath,"report_for_top3.tsv"))
    top3_raw = top3ByGenerator(generator)
    return triq_raw, msstats_raw, msq_raw, top3_raw


for mDir in methodsDir:
    print(f"Processing data from {mDir}")
    subd = os.path.join(dataDir, mDir)
    outd = os.path.join(resDir, mDir)
    raws = readMethods(subd)




df.drop(df[df[pCol].str.contains("DECOY_")].index, inplace=True)
negative = df[pCol].str.contains("_HUMAN").astype(int)
df[negCol] = negative.cumsum()
df[posCol] = (1-negative).cumsum()
df[eCol] = df[negCol]/(df[posCol]+df[negCol])
methodNames = ["Triqler", "MSStats", "MSqRobSum","Top3"]
counted = [countProteins(raw, "Protein", method) for raw, method in zip(raws, methodNames)]
pn = pd.concat([df[[posCol,negCol,"method"]] for df in counted], ignore_index=True, sort=False)

# There is something wrong with the calibration plotting method in pqplot
# negative.cumsum() thing does not work as expected?
# There is something wrong with the plotting in pqplot script 
# Just redo the actual vs amount of differential HeLa.

plt.plot(counted[0]["FDR"], counted[0][eCol])
plt.plot(counted[1]["FDR"], counted[1][eCol])
plt.plot(counted[2]["FDR"], counted[2][eCol])
plt.plot(counted[3]["FDR"], counted[3][eCol])


counted[0]
counted[1]
counted[2]

counted[3]








