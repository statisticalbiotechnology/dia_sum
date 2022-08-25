#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:51:34 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import argparse

def adjust_protein_groups(df):
    protein_list =  df["Protein.Ids"].unique()
    pg_list = df[df["Protein.Ids"].str.contains(";")]["Protein.Ids"].unique()
    pg_keep_list = []
    for pg in pg_list:
         for protein in pg.split(";"):
            if protein not in protein_list:
                if pg not in pg_keep_list:
                    pg_keep_list.append(pg)
    df_non_pgs = df[~df["Protein.Ids"].str.contains(";")]
    df_pgs = df[df["Protein.Ids"].str.contains(";")]
    df_pgs = df_pgs[df_pgs["Protein.Ids"].isin(pg_keep_list)]
    return pd.concat([df_non_pgs, df_pgs])

def convert_diann_to_triqler_PXD031322(filename, output, conditions = "ST|LT|Ctrl", fdr_max = 1.00, protein_group_adjustment = True):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "Q.Value",
                                "Precursor.Normalised", "Stripped.Sequence",
                                "Protein.Ids"])
    if protein_group_adjustment == True:
        df = adjust_protein_groups(df)
    df = df[df["Q.Value"] < fdr_max]

    df["condition"] = df.Run.map(lambda x: x.split("-")[-1])
    
    #df = df[df["condition"].str.contains("LT|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|LT|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|LT")] # Change here for group
    df = df[df["condition"].str.contains(conditions)] # Change here for group

    df["condition"] = df.condition.map(lambda x:x[:-1])
    df.rename({"Run":"run", "Precursor.Charge":"charge", "Q.Value":"searchScore",
               "Precursor.Normalised":"intensity", "Stripped.Sequence":"peptide",
               "Protein.Ids":"proteins"}, axis = 1, inplace = True)
    df = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    df["searchScore"] = -np.log(df["searchScore"])
    df = df.dropna()
    df.to_csv(output, sep = "\t", index = False)





if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots a grid plot of histograms.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--filename', type=str, default = "report.tsv")
    parser.add_argument('--fdr_threshold', type=float,
                        help='Fdr threshold.',
                        default=0.05)
    
    parser.add_argument('--output', type=str,
                        help='output name.')
    
    
    # parse arguments from command line
    args = parser.parse_args()

    filename = args.filename
    fdr_threshold = args.fdr_threshold
    output = args.output

    
    convert_diann_to_triqler_PXD031322(filename = filename, 
                             output = output,
                             fdr_max = fdr_threshold)



    print("Done!")



for i in np.arange(0,1, 0.05):
    convert_diann_to_triqler_PXD031322(filename = "report.tsv", 
                             output = f"triqler_input_report_{round(i,2)}.tsv",
                             fdr_max = round(i,2))
    convert_diann_to_triqler_PXD031322(filename = "report.tsv", 
                             output = f"triqler_input_report_{1}.tsv",
                             fdr_max = 1)

convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_LT.tsv",
                                   fdr_max = 1)

convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_LT_ST.tsv",
                                   fdr_max = 1)


#####



convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_ST_LT_CTRL_adj_pg.tsv",
                                   conditions = "ST|LT|Ctrl",
                                   fdr_max = 1)

convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_ST_CTRL.tsv",
                                   conditions = "ST|Ctrl",
                                   fdr_max = 1)

convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_LT_CTRL.tsv",
                                   conditions = "LT|Ctrl",
                                   fdr_max = 1)

convert_diann_to_triqler_PXD031322(filename = "report.tsv",
                                   output = "triqler_input_report_LT_ST.tsv",
                                   conditions = "ST|LT",
                                   fdr_max = 1)

#import os
#os.chdir("/home/ptruong/data/PXD004873_quick_test")

#filename = "E1603180948_matrix.tsv"
#df = pd.read_csv(filename, sep = "\t")
#df[df.Peptide.str.contains("DECOY")]
#df.columns.
# 453187


#df_ = df[df["Q.Value"] < 0.05]
#df_["Protein.Ids"].map(str.split)
#df["Fragment.Quant.Raw"]



#140970/674416

#250/0.209
