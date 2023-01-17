#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:45:45 2022

@author: ptruong
"""


import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from parsers.parse_triqler import  parse_triqler
import argparse
import numpy.random as npr
import numpy as np

sns.set_context("talk")

negCol = "Differential Protein"
posCol = "Differential Entrapment Protein"
eCol = "Actual Error Rate"


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


def remove_decoy(df, protein_column):   
    res = df[~df[protein_column].str.contains("DECOY_")].copy(deep=True)
    return res

def countProteins(df, method):
    df = remove_decoy(df, protein_column = "Protein")
    df.dropna(subset=["FDR"], inplace = True)
    df.sort_values(by = "FDR", inplace = True)
    negative = df["Protein"].str.contains("Random").astype(int).copy(deep=True)
    df[negCol] = negative.cumsum().copy(deep=True)
    df[posCol] = (1-negative).cumsum().copy(deep=True)
    #if specie == "all":
    #    df[posCol] = (1-negative).cumsum().copy(deep=True)
    #elif specie == "ecoli":
    #    df = df[~df["Protein"].str.contains("_YEAST")]
    #    positive = df["Protein"].str.contains("_ECOLI").astype(int).copy(deep=True)
     #   df[posCol] = positive.cumsum().copy(deep=True)
    #elif specie == "yeast":
    #    df = df[~df["Protein"].str.contains("_ECOLI")]
    #    positive = df["Protein"].str.contains("_YEAST").astype(int).copy(deep=True)
    #    df[posCol] = positive.cumsum().copy(deep=True)
    df["method"] = method
    return df


def read_in_files(triqler_file = "triqler_results/fc_0.96",
                  top3_file = "top3_output_diann.csv",
                  msstats_file = "msstat_output.csv",
                  msqrob2_file = "msqrob2_results.tsv"):
    triqler_results = parse_triqler(triqler_file)
    top3_results = pd.read_csv(top3_file, sep = "\t")
    msstats_results = pd.read_csv(msstats_file, sep = ",")
    msqrob2_results = pd.read_csv(msqrob2_file, sep = ",").rename({"Unnamed: 0":"Protein"},axis=1)
    
    #Rename protein, fdr to same name
    triqler_results = triqler_results.rename({"q_value":"FDR", "protein":"Protein", "log2_fold_change": "log2FC"}, axis = 1)
    top3_results = top3_results.rename({"q":"FDR", "ProteinName":"Protein", "log2(A,B)": "log2FC"}, axis = 1)
    msstats_results = msstats_results.rename({"adj.pvalue":"FDR", "pvalue":"p"}, axis = 1)
    msqrob2_results = msqrob2_results.rename({"adjPval":"FDR", "logFC":"log2FC", "pval":"p"}, axis = 1)
    
    methods = ["Triqler", "Top3", "MsStats", "MsqRob2"]
    data = [triqler_results, top3_results, msstats_results, msqrob2_results]
    zipped_files = zip(methods, data)
    return zipped_files


def get_differential_abundance_count(zipped_files, fc_threshold = 0.01, fdr_threshold = 1.00):
    dfs = []
    test = []
    for method, df in zipped_files:
        if method != "Triqler":
            if fc_threshold > 0:
                df = df[abs(df.log2FC) > fc_threshold]
                df["FDR"] = qvalues(df, pcol = "p")["q"].copy()
        df = df[df.FDR < fdr_threshold].copy()
        df_count = countProteins(df, method)
        test.append(df)
        dfs.append(df_count.loc[:,["Protein", "FDR", posCol, negCol, "method"]])
    res = pd.concat(dfs)
    res = res.reset_index().drop("index", axis = 1)
    return res


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots differential HeLa vs differential non-HeLa lineplot.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--triqler_input', type=str,
                        help='Triqler results file.')

    parser.add_argument('--top3_input', type=str,
                        help='Top3 results file.')

    parser.add_argument('--msstats_input', type=str,
                        help='MsStats results file.')

    parser.add_argument('--msqrob2_input', type=str,
                        help='MSqRob2 results file.')

    parser.add_argument('--specie', type=str,
                        help='specie "all", "ecoli" or "yeast" as y-axis for the differential plot.',
                        default = "all")
    
    parser.add_argument('--fc_threshold', type=float, default = 0,
                        help='Apply fold-change threshold to Top3, MSstats and MSqRob2.')

    parser.add_argument('--fdr_threshold', type=float, default = 1.00,
                        help='Apply q-value threshold.')

    parser.add_argument('--output', type=str,
                        help='Output name.')

    # parse arguments from command line
    args = parser.parse_args()
    triqler_file = args.triqler_input
    top3_file = args.top3_input
    msstats_file = args.msstats_input
    msqrob2_file = args.msqrob2_input
    specie = args.specie
    fc_threshold = args.fc_threshold
    fdr_threshold = args.fdr_threshold
    output = args.output

    
    print(f"""
          Reading in files:
              triqler : {triqler_file}
              top3 : {top3_file}
              msstats : {msstats_file}
              msqrob2 : {msqrob2_file}
          """)
    zipped_files = read_in_files(triqler_file = triqler_file,
                      top3_file = top3_file,
                      msstats_file = msstats_file,
                      msqrob2_file = msqrob2_file)
    print("Counting differentially abundant proteins")
    df_count = get_differential_abundance_count(zipped_files, fc_threshold = fc_threshold)
    
    print("Plotting lineplot")
    fig, ax = plt.subplots(1, 1, figsize=(18,12))

    #if specie == "all":
    sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
    ax.set_ylabel("Differential Proteins", fontsize = 42)
    #elif specie == "yeast":
    #    sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
    #    ax.set_ylabel("Differential Yeast", fontsize = 42)
    #elif specie == "ecoli":
    #    sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
    #    ax.set_ylabel("Differential E.coli", fontsize = 42)
    ax.set_xlim(-1,200)

    ax.set_xlabel("Differential Entrapment Proteins", fontsize=42)
    #ax.set_ylabel(fontsize=38)

    ax.tick_params(axis='x', which='major', labelsize=44)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=44)
    ax.set_xlim([0,50])

    plt.setp(ax.get_legend().get_texts(), fontsize='38') # for legend text
    plt.setp(ax.get_legend().get_title(), fontsize='42') # for legend title

    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)
    
    

#for method, df in zipped_files:
#    print(method)
    
    
    
    
    
    
