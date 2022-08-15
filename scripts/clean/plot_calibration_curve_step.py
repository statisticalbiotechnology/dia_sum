#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 14:53:50 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
from parsers.parse_triqler import  parse_triqler
import seaborn as sns
import matplotlib.pyplot as plt
import numpy.random as npr
import argparse
sns.set_context("talk")



def bootstrap(invec):
    idx = npr.randint(0, len(invec), len(invec))
    return [invec[i] for i in idx]

def estimatePi0(p, numBoot=100, numLambda=100, maxLambda=0.95):
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

def qvalues(pvalues, pcol = "p"):
    m = pvalues.shape[0] # The number of p-values
    pvalues.sort_values(pcol,inplace=True) # sort the pvalues in acending order
    pi0 = estimatePi0(list(pvalues[pcol].values))
    print("pi_0 estimated to " + str(pi0))
    
    # calculate a FDR(t) as in Storey & Tibshirani
    num_p = 0.0
    for ix in pvalues.index:
        num_p += 1.0
        fdr = pi0*pvalues.loc[ix,pcol]*m/num_p
        pvalues.loc[ix,"q"] = fdr
    
    # calculate a q(p) as the minimal FDR(t)
    old_q=1.0
    for ix in reversed(list(pvalues.index)):
        q = min(old_q,pvalues.loc[ix,"q"])
        old_q = q
        pvalues.loc[ix,"q"] = q
    return pvalues


def pq_data(df, fc):
    df_gt = df[df.log2FC > fc]
    df_lt = df[df.log2FC < -fc]
    df_fc = pd.concat([df_gt, df_lt])
    if "p" in df.columns:
        df_fc["FDR"] = qvalues(df_fc, pcol = "p") ["q"]
    n_array = []
    for q in np.arange(0,0.101,0.001):
        n = (df_fc["FDR"] <= q).sum()
        n_array.append(n)
    res = pd.DataFrame(n_array, index = np.arange(0, 0.101, 0.001), columns = ["DE"])
    return res

def get_pq_data_triqler(triqler, q_val):
    fc_tresh = []
    n_hs = []
    n_ye = []
    n_ec = []
    ecoli_scaling_factor = []
    yeast_scaling_factor = []
    human_scaling_factor = []
   
    #fc = float(file.split("_")[1])
    triqler["specie"] = triqler.protein.map(lambda x:x.split("_")[1])
    df_triq = triqler
    ecoli_factor = (df_triq["specie"] == "ECOLI").sum()/(df_triq["specie"] == "HUMAN").sum()
    yeast_factor = (df_triq["specie"] == "YEAST").sum()/(df_triq["specie"] == "HUMAN").sum()
    human_factor = (df_triq["specie"] == "HUMAN").sum()/(df_triq["specie"] == "HUMAN").sum()
    df_triq = df_triq[df_triq.q_value < q_val]
    n_hs.append((df_triq["specie"] == "HUMAN").sum())
    n_ye.append((df_triq["specie"] == "YEAST").sum())
    n_ec.append((df_triq["specie"] == "ECOLI").sum())
    ecoli_scaling_factor.append(ecoli_factor) #Should be same for all values because we want the full length of protein list
    yeast_scaling_factor.append(yeast_factor) #Should be same for all values because we want the full length of protein list
    human_scaling_factor.append(human_factor)
    #fc_tresh.append(fc)
    df =pd.DataFrame(np.array([n_hs, n_ye, n_ec,
                               human_scaling_factor, yeast_scaling_factor, ecoli_scaling_factor]).T, 
                     columns = ["HUMAN", "YEAST", "ECOLI", "HUMAN_factor", "YEAST_factor", "ECOLI_factor"])
    return df

def get_DE_for_triqler(triqler):
    qs = []
    dfs = []
    q_range = np.arange(0,0.101, 0.001) 
    for q in q_range:
        qs.append(q)
        dfs.append(get_pq_data_triqler(triqler, q))
    
    df_res = pd.concat(dfs).set_index(q_range)
    #df_res = pd.DataFrame(vals, index = np.arange(0,0.101, 0.001), columns = ["HUMAN", "YEAS8", "ECOLI", "HUMAN_factor", "YEAST_factor", "ECOLI_factor"])
    return df_res

def calibration_plot(df, output, xlim = [0,0.10], ylim = [0,0.20]):
    fig, axs = plt.subplots(1, 1, figsize=(10,6))
    #sns.lineplot(x = "FDR", y = "actual_error", data = df, ax = axs, hue = "method")
    sns.lineplot(x = "FDR", y = "Fraction_HeLa", data = df, ax = axs, hue = "Method")
    
    #axs.plot(df[df["Method"] == "Triqler"].FDR, df[df["Method"] == "Triqler"].Fraction_HeLa)
    #axs.plot(df[df["Method"] == "Top3"].FDR, df[df["Method"] == "Top3"].Fraction_HeLa)
    #axs.plot(df[df["Method"] == "MsStats"].FDR, df[df["Method"] == "MsStats"].Fraction_HeLa)
    #axs.plot(df[df["Method"] == "MsqRob2"].FDR, df[df["Method"] == "MsqRob2"].Fraction_HeLa)     

    #axs.set_xlabel(r"\textit{q}-value / FDR ", fontsize=34)
    #axs.set_ylabel("Fraction HeLa", fontsize=38)

    axs.set_xlabel("Estimated FDR/q-value", fontsize=24)
    axs.set_ylabel(r"Fraction of HeLa samples", fontsize=24)

    axs.tick_params(axis='x', which='major', labelsize=21)#labelrotation=90)
    axs.tick_params(axis='y', which='major', labelsize=21)
    axs.set_xlim(xlim)
    axs.set_ylim(ylim)

    def abline(slope, intercept):
        """Plot a line from slope and intercept"""
        #axes = plt.gca()
        x_vals = np.array(axs.get_xlim())
        y_vals = intercept + slope * x_vals
        axs.plot(x_vals, y_vals, 'k--', alpha = 0.7)
    abline(1,0)
    print("Outputting : " + output)
    plt.savefig(output, bbox_inches="tight")

def get_fraction_hela(df_in, method):
    df=df_in.copy()
    df.sort_values(by = "FDR", inplace = True)
    df["count_HUMAN"] = df.Protein.str.contains("_HUMAN").astype(int)
    df["count_ECOLI"] = df.Protein.str.contains("_ECOLI").astype(int)
    df["count_YEAST"] = df.Protein.str.contains("_YEAST").astype(int)
    df["cumsum_HUMAN"] = df["count_HUMAN"].cumsum()
    df["cumsum_ECOLI"] = df["count_ECOLI"].cumsum()
    df["cumsum_YEAST"] = df["count_YEAST"].cumsum()
    df["Fraction_HeLa"] = df["cumsum_HUMAN"] / (df["cumsum_HUMAN"] + df["cumsum_ECOLI"] + df["cumsum_YEAST"])
    df["Method"] = method
    df.fillna(1, inplace = True) # assume NaN is maxed out FDR
    #df = df.reset_index().drop("index",axis = 1)
    return df[["FDR", "Fraction_HeLa", "Method"]]

def get_fraction_hela_df(triqler, top3, msstats, msqrob2):
    fdrs = []
    for df, method in zip([triqler, top3, msstats, msqrob2], ["Triqler", "Top3", "MSstats", "MSqRob2"]):
        fdrs.append(get_fraction_hela(df, method))
    df = pd.concat(fdrs)
    df = df.reset_index().drop("index", axis = 1).copy()
    return df

def threshold_fc(df, fc_threshold):
    fc_threshold = 0.0
    df_fc = df[abs(df["log2FC"]) > fc_threshold].copy()
    df_fc["FDR"] = qvalues(df_fc, pcol = "p")["q"] #recompute FDR for thresholded set
    return df_fc


def main(triqler_file, top3_file, msstats_file, msqrob2_file, output, fc_threshold = 0, xlim = [0,0.10], ylim = [0,0.2]):
    
    #fc_threshold = 0.0 #should be same as triqler fc_eval

    #triqler_file = "triqler_results.csv"
    #top3_file = "top3_results.csv"
    #msstats_file = "msstats_results.csv"
    #msqrob2_file = "msqrob2_results.csv"
    
    triqler = parse_triqler(triqler_file)
    #triqler[~triqler.protein.str.contains("DECOY")]
    triqler["specie"] = triqler.protein.map(lambda x:x.split("_")[1])
    triqler["FDR"] = triqler["q_value"]
    triqler.rename({"protein":"Protein"}, axis = 1, inplace = True)
    top3 = pd.read_csv(top3_file, sep = "\t").rename({"q":"FDR", 'log2(A,B)':"log2FC"}, axis = 1)
    top3.rename({"ProteinName":"Protein"}, axis = 1, inplace = True)
    msstats = pd.read_csv(msstats_file, sep = ",").rename({"adj.pvalue":"FDR", "pvalue":"p"}, axis = 1)
    msstats["specie"] = msstats.Protein.map(lambda x:x.split("_")[1])
    msqrob2 = pd.read_csv(msqrob2_file, sep = ",").rename({"Unnamed: 0":"protein", "adjPval":"FDR", "pval":"p", "logFC":"log2FC"}, axis = 1)
    msqrob2.rename({"protein":"Protein"}, axis = 1, inplace = True)
    msqrob2["specie"] = msqrob2.Protein.map(lambda x:x.split("_")[1])
    msqrob2.rename({"protein":"Protein"}, axis = 1, inplace = True)

    if fc_threshold != 0:
        #triqler = threshold_fc(triqler, fc_threshold)
        top3 = threshold_fc(top3, fc_threshold)
        msstats = threshold_fc(msstats, fc_threshold)
        msqrob2 = threshold_fc(msqrob2, fc_threshold)
        
    df = get_fraction_hela_df(triqler, top3, msstats, msqrob2)
    
    calibration_plot(df = df, output = output, xlim = xlim, ylim = ylim)


                            

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots calibration plots.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--triqler_input', type=str,
                        help='Triqler results file.')

    parser.add_argument('--top3_input', type=str,
                        help='Top3 results file.')

    parser.add_argument('--msstats_input', type=str,
                        help='MsStats results file.')

    parser.add_argument('--msqrob2_input', type=str,
                        help='MSqRob2 results file.')

    parser.add_argument('--fc_threshold', type=float,
                        help='Log2FC threshold to apply on Top3, MsStats and MSqRob2 results. NOTE: Triqler does not need Log2FC thresholding. Use the same fold_change_eval parameter when running triqler.')

    parser.add_argument('--output', type=str,
                        help='Output name.')

    parser.add_argument('--x_min', type=float, default = 0.0,
                        help='x_min for plot.')

    parser.add_argument('--x_max', type=float, default = 0.1,
                        help='x_max for plot.')

    parser.add_argument('--y_min', type=float, default = 0.0,
                        help='y_min for plot.')

    parser.add_argument('--y_max', type=float, default = 0.2,
                        help='y_max for plot.')


    # parse arguments from command line
    args = parser.parse_args()
    triqler_file = args.triqler_input
    top3_file = args.top3_input
    msstats_file = args.msstats_input
    msqrob2_file = args.msqrob2_input
    fc_threshold = args.fc_threshold
    output = args.output
    x_min = args.x_min
    x_max = args.x_max
    y_min = args.y_min
    y_max = args.y_max
    
    main(triqler_file = triqler_file, 
         top3_file = top3_file, 
         msstats_file = msstats_file, 
         msqrob2_file = msqrob2_file, 
         output = output, 
         fc_threshold = fc_threshold, 
         xlim = [x_min, x_max], 
         ylim = [y_min, y_max])