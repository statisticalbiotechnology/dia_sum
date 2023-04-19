#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:19:56 2022

@author: ptruong
"""
import pandas as pd 
import os 

os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/supplementary_material")

ST_ctrl = pd.read_excel("S2_list_of_proteins_differentially_regulated_between_ST_and_Ctrl.xlsx", header = 1)
LT_ctrl = pd.read_excel("S3_list_of_proteins_differentially_regulated_between_LT_and_Ctrl.xlsx", header = 1)
diffReg = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
keggEnrichment = pd.read_excel("S5_KEGG_enrichment_analysis_of_proteins_in_six_subclusters.xlsx", header = 1)
gsvaEnrichment = pd.read_excel("S6_GSVA_enrichment_analysis.xlsx", header = 1)

# Parse their data and make the clusters thier way.

ST_ctrl[~ST_ctrl["Condition"].isin(["Not_sig"]) & (ST_ctrl['log2（FC）'] > 1.5) ] #12
ST_ctrl[~ST_ctrl["Condition"].isin(["Not_sig"]) & (ST_ctrl['log2（FC）'] > 0.67)] #70
12+70

pd.concat([ST_ctrl[~ST_ctrl["Condition"].isin(["Not_sig"]) & (ST_ctrl['log2（FC）'] > 1.5) ], 
ST_ctrl[~ST_ctrl["Condition"].isin(["Not_sig"]) & (ST_ctrl['log2（FC）'] > 0.67)]])

LT_ctrl[~LT_ctrl["Condition"].isin(["Not_sig"]) & (LT_ctrl['log2（FC）'] > 1.5) ] #61
LT_ctrl[~LT_ctrl["Condition"].isin(["Not_sig"]) & (LT_ctrl['log2（FC）'] > 0.67)] #251
61+251

pd.concat([LT_ctrl[~LT_ctrl["Condition"].isin(["Not_sig"]) & (LT_ctrl['log2（FC）'] > 1.5) ],
LT_ctrl[~LT_ctrl["Condition"].isin(["Not_sig"]) & (LT_ctrl['log2（FC）'] > 0.67)]]) #251

#280
diffReg.columns
diffReg = diffReg[~diffReg["Condition"].isin(["Others"])]

lt_st = pd.concat([diffReg[(diffReg["LT-ST_diff"] > 1.5) & (diffReg["LT-ST_p adj"] < 0.05)], diffReg[(diffReg["LT-ST_diff"] < 0.67) & (diffReg["LT-ST_p adj"] < 0.05)] ])
st_ctrl = pd.concat([diffReg[(diffReg["ST-Ctrl_diff"] > 1.5) & (diffReg["ST-Ctrl_p adj"] < 0.05)], diffReg[(diffReg["ST-Ctrl_diff"] < 0.67) & (diffReg["ST-Ctrl_p adj"] < 0.05)] ])
lt_ctrl = pd.concat([diffReg[(diffReg["LT-Ctrl_diff"] > 1.5) & (diffReg["LT-Ctrl_p adj"] < 0.05)], diffReg[(diffReg["LT-Ctrl_diff"] < 0.67) & (diffReg["LT-Ctrl_p adj"] < 0.05)] ])


st_ctrl[st_ctrl["Condition"].isin(["C5"])] # 114
st_ctrl[st_ctrl["Condition"].isin(["C6"])] # 244
st_ctrl[st_ctrl["Condition"].isin(["C5", "C6"])] #358

lt_st[lt_st["Condition"].isin(["C3"])] #137
lt_st[lt_st["Condition"].isin(["C4"])] #422
lt_st[lt_st["Condition"].isin(["C3", "C4"])] #559



lt_st[lt_st["Condition"].isin(["C3", "C4"])] #559
st_ctrl[st_ctrl["Condition"].isin(["C5","C6"])] #358
lt_ctrl[lt_ctrl["Condition"].isin(["C1", "C2"])] #87
















