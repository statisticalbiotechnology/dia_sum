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

sns.set_context("talk")

negCol = "Differential HeLa"
posCol = "Differential non-HeLa"
eCol = "Actual Error Rate"

def remove_decoy(df, protein_column):   
    res = df[~df[protein_column].str.contains("DECOY_")].copy(deep=True)
    return res

def countProteins(df, method, specie = "all"):
    df = remove_decoy(df, protein_column = "Protein")
    df.dropna(subset=["FDR"], inplace = True)
    df.sort_values(by = "FDR", inplace = True)
    negative = df["Protein"].str.contains("_HUMAN").astype(int).copy(deep=True)
    df[negCol] = negative.cumsum().copy(deep=True)
    if specie == "all":
        df[posCol] = (1-negative).cumsum().copy(deep=True)
    elif specie == "ecoli":
        df = df[~df["Protein"].str.contains("_YEAST")]
        positive = df["Protein"].str.contains("_ECOLI").astype(int).copy(deep=True)
        df[posCol] = positive.cumsum().copy(deep=True)
    elif specie == "yeast":
        df = df[~df["Protein"].str.contains("_ECOLI")]
        positive = df["Protein"].str.contains("_YEAST").astype(int).copy(deep=True)
        df[posCol] = positive.cumsum().copy(deep=True)
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
    triqler_results = triqler_results.rename({"q_value":"FDR", "protein":"Protein"}, axis = 1)
    top3_results = top3_results.rename({"q":"FDR", "ProteinName":"Protein"}, axis = 1)
    msstats_results = msstats_results.rename({"adj.pvalue":"FDR"}, axis = 1)
    msqrob2_results = msqrob2_results.rename({"adjPval":"FDR"}, axis = 1)
    
    methods = ["Triqler", "Top3", "MsStats", "MsqRob2"]
    data = [triqler_results, top3_results, msstats_results, msqrob2_results]
    zipped_files = zip(methods, data)
    return zipped_files


def get_differential_abundance_count(zipped_files, specie = "all"):
    dfs = []
    for method, df in zipped_files:
        df_count = countProteins(df, method, specie = specie)
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

    parser.add_argument('--output', type=str,
                        help='Output name.')

    # parse arguments from command line
    args = parser.parse_args()
    triqler_file = args.triqler_input
    top3_file = args.top3_input
    msstats_file = args.msstats_input
    msqrob2_file = args.msqrob2_input
    specie = args.specie
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
    df_count = get_differential_abundance_count(zipped_files, specie = specie)
    
    print("Plotting lineplot")
    fig, ax = plt.subplots(1, 1, figsize=(18,12))

    if specie == "all":
        sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
        ax.set_ylabel("Differential non-HeLa", fontsize = 38)
    elif specie == "yeast":
        sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
        ax.set_ylabel("Differential yeast", fontsize = 38)
    elif specie == "ecoli":
        sns.lineplot(x = negCol, y = posCol, hue = "method", data = df_count, ax = ax, linewidth = 5)
        ax.set_ylabel("Differential ecoli", fontsize = 38)
    ax.set_xlim(-1,200)

    ax.set_xlabel("Differential HeLa", fontsize=34)
    #ax.set_ylabel(fontsize=38)

    ax.tick_params(axis='x', which='major', labelsize=42)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=42)
    #ax.set_xlim(xlim)

    plt.setp(ax.get_legend().get_texts(), fontsize='32') # for legend text
    plt.setp(ax.get_legend().get_title(), fontsize='38') # for legend title

    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)
    
    

#for method, df in zipped_files:
#    print(method)
    
    
    
    
    
    