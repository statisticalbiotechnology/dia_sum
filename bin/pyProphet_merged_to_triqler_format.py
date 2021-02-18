#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 17:05:05 2021

@author: ptruong
"""

import pandas as pd
 

df = pd.read_csv("merged_recolumned.tsv", sep = "\t")

run = df.run_id
condition = df.filename
charge = df.Charge
intensity = df.Intensity
searchScore = df.m_score_peptide_run_specific
peptide = df.FullPeptideName
protein = df.ProteinName


df_triq = df[["run_id", "filename", "Charge", "Intensity", "m_score_peptide_run_specific", "FullPeptideName", "ProteinName"]]

df_triq = df_triq.rename(columns={"run_id":"run", 
                        "filename":"condition",
                        "Charge":"charge", 
                        "Intensity":"intensity", 
                        "m_score_peptide_run_specific":"searchScore", 
                        "FullPeptideName":"peptide", 
                        "ProteinName":"proteins"})
df_triq.to_csv("triqler_formatted.csv", sep = "\t", index=False)





#df_triq[df_triq.run == "BR_Repl2_TR_R01_10"].condition.unique()


#for i in df_triq.run.unique():
#    print(len(df_triq[df_triq.run == "BR_Repl2_TR_R01_10"].condition.unique()))

