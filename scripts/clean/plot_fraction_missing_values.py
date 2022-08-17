#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 11:25:57 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
import argparse


sns.set_context("talk")

#import warnings
#warnings.filterwarnings("ignore")

def read_file(triqler_input_file):
    df = pd.read_csv(triqler_input_file, sep = "\t")
    df = df[df['proteins'].notna()]
    decoy_mapper = lambda x: x.split("_")[0]
    df["decoy"] = (df.proteins.map(decoy_mapper) == "DECOY").copy(deep=True)
    df = df[df["decoy"] == False] #Filter away decoy peptides and protein
    return df

def get_peptide_intensity_count(df):
    return df.groupby(["condition", "peptide"]).count().intensity 

def split_peptide_intensity_count_into_conditions(df_peptide_intensity_count, condition):
    return df_peptide_intensity_count[(df_peptide_intensity_count<3) 
                                  & (df_peptide_intensity_count.index.get_level_values("condition").isin([condition]))].index.get_level_values("peptide")

def create_pivot_table_with_missing_count_per_condition(df):
    print("Computing peptide intensity count...")
    df_peptide_intensity_count = get_peptide_intensity_count(df)
    pivot_missing_list = []
    conditions_cols = {}
    for i in df_peptide_intensity_count.index.get_level_values("condition").unique():
        # We split up between conditions because we do not want samples to partially exist in condition 1 and 2 
        # and count them as full samples.
        # We check which peptides has lower than count 3 i.e. one for each sample for each condition>
        cond_missing_value_peptides = split_peptide_intensity_count_into_conditions(df_peptide_intensity_count, condition = i)
        df_condition_missing = df[df.peptide.isin(cond_missing_value_peptides) & (df.condition.isin([i]))]
        # We create a pivot table (rows = peptides, columns = samples)
        pivot_condition_missing = pd.pivot_table(df_condition_missing, values = "intensity", index = ["condition", "peptide"], columns = "run")
        pivot_missing_list.append(pivot_condition_missing)
        conditions_cols[i] = pivot_condition_missing.columns
        
    pivot_missing = pd.concat(pivot_missing_list, axis = 0)
    
    # Compute mean, standard deviation and na_count
    pivot_missing["mu"] = pivot_missing.mean(axis = 1)
    pivot_missing["std"] = pivot_missing.std(axis = 1) 
    pivot_missing["na_count"] = pivot_missing.isna().sum(axis = 1) - 3 # remove three because we do not want to count missing on both samples
    return pivot_missing, conditions_cols


def get_condition_mean_imputed_nans(pivot_missing, conditions_cols):
    df_cond_imputed_mean_nans_list = []
    for i in pivot_missing.index.get_level_values("condition").unique():
        # We fill na and create nan_bool, so we can hide the filled in values again, and a vals_bool so we can we 
        # the non-imputed values.
        df_cond = pivot_missing[pivot_missing.index.get_level_values("condition") == i][list(conditions_cols[i])]
        vals_bool = ~df_cond.isna()
        #vals_bools = vals_bool.replace(False, np.nan)
        nan_bool = df_cond.isna()
        nan_bool = nan_bool.replace(False, np.nan)
        df_cond = df_cond.T.fillna(df_cond.mean(axis=1)).T #row average fillna
    
        # Keep imputed nans
        df_cond_imputed_mean_nans = (df_cond*nan_bool).droplevel(level = "condition").melt().dropna()
        # Keep values 
        #df_cond_vals = (df_cond*vals_bool).droplevel(level = "condition").melt().dropna()
        df_cond_imputed_mean_nans_list.append(df_cond_imputed_mean_nans)
     
    #Create array for imputed means. 
    
    df_imputed_nans_with_condition_mean = pd.concat(df_cond_imputed_mean_nans_list).reset_index().drop("index", axis = 1)
    return df_imputed_nans_with_condition_mean

def compute_missing_value_factions(df, bins = np.arange(0,1000,10)):
    #bins = np.arange(0,1000,10)
    #df_peptide_intensity_count = get_peptide_intensity_count(df)
    
    # Create array for intensities.
    df_intensities = df.intensity

    pivot_missing, conditions_cols = create_pivot_table_with_missing_count_per_condition(df)

    df_imputed_nans_with_condition_mean = get_condition_mean_imputed_nans(pivot_missing, conditions_cols)
    
    # Bin the imputed values and intensities.
    df_binned_imputed_mean_nans = pd.cut(df_imputed_nans_with_condition_mean["value"], bins, include_lowest=True)
    df_binned_vals = pd.cut(df_intensities, bins, include_lowest=True)
    
    # For each bin  count(imputed) / count(intensities)
    # Now each bin is an intensity range and each ratio is the fraction of missing value for that intensity.
    df_binned_missing_value_fraction = (df_binned_imputed_mean_nans.value_counts() / df_binned_vals.value_counts())
    df_binned_missing_value_fraction = pd.DataFrame(df_binned_missing_value_fraction.values, index = df_binned_missing_value_fraction.index, columns = ["fraction"]).reset_index()
    df_binned_missing_value_fraction.index = np.arange(0, 990, 10)
    return df_binned_missing_value_fraction

# From triqler code
# logit is from hyperparameters.py
def logit(x, muLogit, sigmaLogit):
    return 0.5 + 0.5 * np.tanh((np.array(x) - muLogit) / sigmaLogit)

# pMissing is from pgm.py
def pmissing(x, muLogit, sigmaLogit):
    return 1-logit(x, muLogit, sigmaLogit) 

parser = argparse.ArgumentParser(
    description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots differential HeLa vs differential non-HeLa lineplot.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--triqler_input', type=str,
                    help='Triqler input file (NOTE: Not the results file).')


parser.add_argument('--output', type=str,
                    help='Output file name.')


# parse arguments from command line
args = parser.parse_args()
triqler_input_file = args.triqler_input
output = args.output

def main():
    print(f"Reading in: {triqler_input_file}")
    df = read_file(triqler_input_file)
    print("Imputing nans with mean...")

    print("Computing missing value fraction...")
    df_binned_missing_value_fraction = compute_missing_value_factions(df, bins = np.arange(0,1000,10))
    
    print("Generating plot...")
    fig, ax = plt.subplots(1, 1, figsize=(16,12))
    xdata = df_binned_missing_value_fraction.index
    ydata = df_binned_missing_value_fraction.fraction
    ax.plot(xdata, ydata, 'b-', label='Fraction missing values', linewidth = 5)

    # We fit the fraction data we have to pmissings
    popt, pcov = curve_fit(pmissing, xdata, ydata)
    ax.plot(xdata, pmissing(xdata, popt[0], popt[1]), "r--", label='fit: muLogit=%5.3f, sigmaLogit=%5.3f' % tuple(popt), linewidth=5)
    
    # NOTE THIS IS NOT LOG-INTENSITY
    ax.set_xlabel("Peptide intensity", fontsize = 34)
    #ax.tick_params(axis='both', which='major', labelsize=20)
    ax.tick_params(axis='x', which='major')
    #ax.set_ylabel("Fraction missing values within binned interval", fontsize = 38)
    ax.set_ylabel("Fraction missing values", fontsize = 34)
    ax.set_xlim(0, 100)
    #ax.set_ylim(0, 0.0)
    ax.legend(fontsize=24)
    ax.tick_params(axis='x', which='major', labelsize=21)#labelrotation=90)
    ax.tick_params(axis='y', which='major', labelsize=21)

    #ax.set_title("DIANN - Fraction Missing Values for mean intensity", fontsize = 22, fontweight = "bold")
    fig = ax.get_figure()
    print(f"Saving output {output}")
    fig.savefig(output)
    print("Done!")    
    
if __name__ == "__main__":
    main()
