#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 09:22:40 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy import stats
import argparse

import os 

from q_value import qvalues

def read_in_diann(file = "report_for_top3.tsv", q_value_threshold = 0.01, log = True):
    df = pd.read_csv(file, sep = "\t")                      
    df = df.rename({"Stripped.Sequence":"peptide", "Precursor.Quantity":"intensity"}, axis = 1)
    if log == True:
        df.intensity = np.log(df.intensity)
    df = df[df["Q.Value"] < q_value_threshold]
    peptides_included_in_all_samples = (df.groupby("peptide").count() >= 6).index
    df = df[df["peptide"].isin(peptides_included_in_all_samples)]
    df = df[df["Protein.Ids"].notna()] # drop proteins not identified
    df = df[~df["Protein.Ids"].str.contains("DECOY")] # drop decoy
    return df

def get_peptide_mu_sigma(df, quantiles = 8):
    df_means = df.groupby("peptide").mean()
    df_stat = pd.DataFrame(df_means.intensity.values, index = df_means.index, columns = ["mu"])
    df_stat["std"] = df.groupby("peptide").std().intensity
    df_stat["cv"] = df_stat["std"] / df_stat["mu"]
    df_stat["quantile_bin_mu"] = pd.qcut(df_stat["mu"], q=quantiles)

    bin_median_function = lambda x: round((x.left + x.right)/2, 2)
    df_stat["quantile_bin_mu_median_of_bin_range"] = df_stat["quantile_bin_mu"].apply(bin_median_function)
    return df_stat

def plot_boxplot(df, log_labels = True, ylim = False, hline = "median", output = "output.png"):
    f, ax = plt.subplots(1, 1, figsize = (12,12))
    #sns.violinplot(x='quantile_bin_mu', y='std', data=df_stats, ax=ax)
    sns.boxplot(x='quantile_bin_mu_median_of_bin_range', y='std', data=df, ax=ax, color = "lightblue")

    if log_labels == True:
        ax.set_ylabel("Standard Deviation of log of Peptide intensity")
        ax.set_xlabel("Log of Peptide")
    else:
        ax.set_ylabel("Standard Deviation of Peptide intensity")
        ax.set_xlabel("Peptide")
    ax.tick_params(axis='x', which='major',labelrotation=90)
    if ylim != False:
        ax.set_ylim(ylim)
    if hline == "mean":
        ax.axhline(df["std"].mean(), ls='--', color = "red")
    elif hline == "median":
        ax.axhline(df["std"].median(), ls='--', color = "red")
    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig("output")

def read_in_osw(file_directory, m_score_threshold, log = True):
    dfs = []
    for i in os.listdir(file_directory):
        df = pd.read_csv(file_directory + i, sep = "\t")
        df = df[df.m_score < m_score_threshold]
        df = df[df["ProteinName"].notna()] # drop proteins not identified
        df = df[df["ProteinName"].str.contains("DECOY")]
        df = df.rename({"Sequence":"peptide", "Intensity":"intensity"}, axis = 1)
        if log == True:
            df.intensity = np.log(df.intensity)
        dfs.append(df)
    df = pd.concat(dfs)
    return df

def main_osw(
    file_directory = '/hdd_14T/data/PXD002952/20210614_dataset/result_files_20220214/DID/osw_results/',
    m_score_threshold = 0.01,
    log = True,
    quantiles = 8,
    ylim = [0,5],
    hline_type = "median",
    output = "output.png"
    ):
    df = read_in_osw(file_directory, m_score_threshold = m_score_threshold, log = log)
    df_stats = get_peptide_mu_sigma(df, quantiles = quantiles).dropna()
    plot_boxplot(df_stats, log_labels = log, ylim = ylim, hline = hline_type, output = output)

def main_diann(
    file_name = "report_for_top3.tsv",
    q_value_threshold = 0.01,
    log = True,
    quantiles = 8,
    ylim = [0,5],
    hline_type = "median",
    output = "output.png"
    ):
    
    df  = read_in_diann(file = file_name, q_value_threshold = 0.01, log = log)
    df_stats = get_peptide_mu_sigma(df, quantiles = quantiles).dropna()
    plot_boxplot(df_stats, log_labels = log, ylim = ylim, hline = hline_type, output = output)
    
def main(file_directory, threshold, log, quantiles, ylim, hline_type, input_type, output):
    input_type = input_type.upper()
    if input_type == "OSW":
        main_osw(
                file_directory = file_directory,
                m_score_threshold = threshold,
                log = log,
                quantiles = quantiles,
                ylim = ylim,
                hline_type = hline_type,
                output
                )
    elif input_type == "DIANN":
        main_diann(
                file_name = file_directory,
                q_value_threshold = threshold,
                log = log,
                quantiles = quantiles,
                ylim = ylim,
                hline_type = hline_type,
                output
                )
    else:
        print("No input_type specified. Exiting!")
    
    
parser = argparse.ArgumentParser(
    description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots differential HeLa vs differential non-HeLa lineplot.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--file_directory', type=str,
                    help='file_directory should be a directory for "OSW" input_type containing the *dscore.csv output from pyprophet, or the report.tsv output from DIA-NN.')

parser.add_argument('--threshold', type=float,
                    help='threshold is the q_value threshold for "OSW" input_type, or the m_score threshold for "DIANN" input_type.',
                    default=0.01)

parser.add_argument('--log', type=bool,
                    help='True/False, if log == True, log-transform is applied to the intensities.',
                    default=True)

parser.add_argument('--quantiles', type=int,
                    help='<int> number of quantiles for the quantile binning.',
                    default=8)

#CHECK THIS
parser.add_argument('--ylim', type=list,
                    help='y limits for the boxplot.',
                    default=False)

parser.add_argument('--hline_type', type=str,
                    help='hline_type False, "median" or "mean" hline for the boxplot.',
                    default="median")

parser.add_argument('--input_type', type=str,
                    help='"OSW" or "DIANN" input type.')

parser.add_argument('--output', type=str,
                    help='Output name.')

# parse arguments from command line
args = parser.parse_args()
file_directory = args.file_directory
threshold = args.threshold
log = args.log
quantiles = args.quantiles
ylim = args.ylim
hline_type = args.input_type
input_type = args.output

output = args.output

if __name__ == "__main__":
    main(file_directory, threshold, log, quantiles, ylim, hline_type, input_type, output)
    
    
    
    
