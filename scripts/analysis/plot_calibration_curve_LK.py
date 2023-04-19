#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:25:18 2022

@author: ptruong
"""
#!/usr/bin/env python3

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

dataDir = "results"
#resDir = "result_plot"
methodsDir = ["ID", "PS"]
eCol = "Actual Error Rate"
negCol = "Differential HeLa"
posCol = "Differential non-HeLa"


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
    names = ["q_value","posterior_error_prob","protein"] + [str(a) for a in range(80)]

    triq_raw = pd.read_csv(os.path.join(resPath, "triqler_results.csv"),
                                usecols=[0,2], sep='\t', skiprows=1, names=names) 
    triq_raw.drop(triq_raw[triq_raw["protein"].str.contains("DECOY_")].index, inplace=True)
    triq_raw.rename(columns={'q_value':"FDR", "protein":"Protein"}, inplace=True)
    
    top3_raw = pd.read_csv(os.path.join(resPath, "top3_results.csv"), sep = "\t")
    top3_raw.rename(columns={'q':"FDR", "ProteinName":"Protein"}, inplace = True)
    top3_raw.sort_values("p", inplace=True)

    msstats_raw = pd.read_csv(os.path.join(resPath, "msstats_results.csv"),
                                usecols=["Protein","pvalue","adj.pvalue"]) 
    msstats_raw.sort_values("pvalue", inplace=True)
    msstats_raw.rename(columns={'adj.pvalue' : "FDR"}, inplace=True)
    msq_raw = pd.read_csv(os.path.join(resPath, "msqrob2_results.csv"),
                                names=["Protein","logFC","se","df","t","pval","adjPval"], skiprows=1, usecols=["Protein","pval","adjPval"]) 
    msq_raw.dropna(inplace=True)                            
    msq_raw.sort_values("pval", inplace=True)
    msq_raw.rename(columns={'adjPval' : "FDR", "" : "Protein"}, inplace=True)

    #if "ID" in resPath:
    #    generator = oswGenerator(os.path.join(resPath,"osw_results"))
    #else:
     #   generator = specialGenerator(os.path.join(resPath,"report_for_top3.tsv"))
    #top3_raw = top3ByGenerator(generator)
    return triq_raw, top3_raw, msstats_raw, msq_raw, 

def plot_calibration_LK(mDir = "ID"):
    #for mDir in methodsDir:
    print(f"Processing data from {mDir}")
    subd = os.path.join(dataDir, mDir)
    #outd = os.path.join(resDir, mDir)
    raws = readMethods(subd)
    methodNames = ["Triqler", "Top3", "MSStats", "MSqRob2"]
    counted = [countProteins(raw, "Protein", method) for raw, method in zip(raws, methodNames)]
    #print(counted[3].to_string())
    pn = pd.concat([df[[posCol,negCol,"method"]] for df in counted], ignore_index=True, sort=False)
    pn.drop(pn[pn[negCol]>100].index, inplace=True)
    er_fdr = pd.concat([df[[eCol, "FDR", "method"]] for df in counted], ignore_index=True, sort=False)
    er_fdr.drop(er_fdr[er_fdr[eCol]>0.2].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr["FDR"].diff(periods=-1)>0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    er_fdr.drop(er_fdr[(er_fdr[eCol].diff()<0.) & (er_fdr["method"].shift(periods=-1) == er_fdr["method"])].index, inplace=True)
    pn.drop_duplicates(subset=[posCol, 'method'], keep='last', inplace=True)
    pn.drop_duplicates(subset=[negCol, 'method'], keep='last', inplace=True)
    # print(pn.to_string())
    # print(er_fdr.to_string())
    #g = sns.lineplot(data=pn, x=negCol, y=posCol, hue="method", ci=None)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles=handles[1:], labels=labels[1:])
    #g.get_legend().set_title(None)
    #plt.xlim(0.,50.)
    #plt.tight_layout()
    #plt.show()
    #plt.savefig(os.path.join(outd, "pos_neg.png"))
    #plt.close()
    fig, g = plt.subplots(1, 1, figsize=(10,6))
    sns.lineplot(data=er_fdr, x="FDR", y=eCol, hue="method", ci=None, ax = g)
    g.get_legend().set_title(None)
    plt.xlim(0.,0.1); plt.ylim(0.,0.3)
    plt.tight_layout()
    def abline(slope, intercept):
        """Plot a line from slope and intercept"""
        #axes = plt.gca()
        x_vals = np.array(g.get_xlim())
        y_vals = intercept + slope * x_vals
        g.plot(x_vals, y_vals, 'k--', alpha = 0.7)
    abline(1,0)
    xlim = [0, 0.10]
    ylim = [0, 0.2]
    g.set_xlim(xlim)
    g.set_ylim(ylim)
    
    plt.show()
    #plt.savefig(os.path.join(outd, "er_fdr.png"))
    #plt.close()
