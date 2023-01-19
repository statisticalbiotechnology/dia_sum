#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 13:32:59 2023

@author: ptruong
"""

import pandas as pd
import numpy as np
import argparse


def convert_diann_to_triqler_PXD031322(filename, output, fdr_max = 1.00):
    df = pd.read_csv(filename, sep = "\t", 
                     usecols = ["Run", "Precursor.Charge", "Q.Value",
                                "Precursor.Normalised", "Stripped.Sequence",
                                "Protein.Ids"])

    df = df[df["Q.Value"] < fdr_max]

    df["condition"] = df.Run.map(lambda x: x.split("_")[2]).map({"7951":"metastatic",
                                                                  "SH4":"metastatic",
                                                                  "HTB69":"metastatic",
                                                                  "SK":"primary",
                                                                  "A375":"primary",
                                                                  "G361":"primary"})



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
