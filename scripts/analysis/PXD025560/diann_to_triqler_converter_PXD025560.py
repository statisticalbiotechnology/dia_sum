#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:51:34 2022

@author: ptruong
"""

import pandas as pd
import numpy as np
import re
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

def convert_diann_to_triqler_PXD031322(filename, metadata_file, output, fdr_max = 1.00, protein_group_adjustment = True):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "Q.Value",
                                "Precursor.Normalised", "Stripped.Sequence",
                                "Protein.Ids"])
    meta = pd.read_csv(metadata_file, sep = "/")
    
    mapper = meta.set_index("SampleID")["Histology"].to_dict()

    df["sample_id"] =  df.Run.map(lambda x: re.findall('\d+',x.split("_")[-1])[0]).astype(int)
    df["sample_class"] = df.sample_id.map(mapper)
    
    #if protein_group_adjustment == True:
    #    df = adjust_protein_groups(df)
    df = df[df["Q.Value"] < fdr_max]

    df["condition"] = df.sample_class
    
    #df = df[df["condition"].str.contains("LT|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|LT|Ctrl")] # Change here for group
    #df = df[df["condition"].str.contains("ST|LT")] # Change here for group
    #df = df[df["condition"].str.contains(conditions)] # Change here for group

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

