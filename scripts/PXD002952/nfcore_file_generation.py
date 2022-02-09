#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:23:43 2021

@author: ptruong
/hdd_14T/data/PXD002952/20210614_dataset
"""


import os 
import pandas as pd
import numpy as np 

def generate_nfcore_input(filepath):
    # location of DIA .mzML
    files = []
    for file in os.listdir(filepath):
        if file[-5:] == ".mzML":
            files.append(file)
            
    df = pd.DataFrame(files, columns = ["Spectra_Filepath"])  
    
    x = df.iloc[0][0]
    
    sample_mapper = lambda x: x.split("_")[5]
    condition_mapper = lambda x: int(x.split("_")[8])
    bioReplicate_mapper = lambda x: int(x.split("_")[-1].split(".")[0][-1])
    
    df = df.reset_index().rename({"index": "Sample"}, axis = "columns")
    #df["Sample"] += 1
    df["BatchID"] = "LFQ"
    df["BatchID"] = 1
    df["MSstats_Condition"] = df.Spectra_Filepath.map(condition_mapper)
    #df["BatchID"] = df["MSstats_Condition"]
    df["MSstats_BioReplicate"] = df.Spectra_Filepath.map(bioReplicate_mapper)
    df[["Sample", "BatchID", "MSstats_Condition", "MSstats_BioReplicate", "Spectra_Filepath"]]
    return df

def generate_nfcore_spectral_library_input(location = ""):
    headers = ["Sample" , "BatchID", "Library_Filepath"]
    data = ["1", "LFQ", location + "easypqp_lib_openswath_specie.tsv"] # from EasyPQP
    #data = ["1", "1", location + "easypqp_lib_openswath_specie.tsv"] # from EasyPQP
    #data = [["1", "1", location + "easypqp_lib_openswath_specie.tsv"],
    #        ["2", "2", location + "easypqp_lib_openswath_specie.tsv"]]
    df = pd.DataFrame(np.array(data).T, index = headers).T
    return df

def generate_nfcore_irts_input(location = ""):
    headers = ["Sample" , "BatchID", "Library_Filepath"]
    #data = ["1", "LFQ", "hroest_DIA_iRT.TraML"]
    data = ["1", "LFQ", location + "irt.tsv"] #from MSFragger
    #data = ["1", "1", location + "hroest_DIA_iRT.TraML"]
    #data = [["1", "1", location + "irt.tsv"],
    #        ["2", "2", location + "irt.tsv"]]
    #data = [["1", "1", location + "hroest_DIA_iRT.TraML"],
    #        ["2", "2", location + "hroest_DIA_iRT.TraML"]]
    df = pd.DataFrame(np.array(data).T, index = headers).T
    return df

df = generate_nfcore_input(os.getcwd())
df.to_csv("nf_sample_input.csv", sep = "\t", index = False)

location = "diaproteomics_20210707/"
#location = ""
lib_df = generate_nfcore_spectral_library_input(location = location)
lib_df.to_csv("nf_input_spectral_library.csv", sep = "\t", index = False)

irts_df = generate_nfcore_irts_input(location = location)
irts_df.to_csv("nf_input_irt_sheet.csv", sep = "\t", index = False)

def generate_input_sheet_dda(mzml_dir = "", pepxml_dir = ""):
    headers = ["Sample", "BatchID", "Spectra_Filepath", "Id_Filepath"]    
    dda_files = []
    pepxml_files = []
    for file in os.listdir(mzml_dir):
        if file[-5:] == ".mzML": #MSFragger output
            dda_files.append(mzml_dir + f"{file}")
    for file in os.listdir(pepxml_dir):
        #if file[-7:] == "pep.xml": #MSFragger output
        if file[-7:] == ".pepXML":
            pepxml_files.append(pepxml_dir + f"{file}")
    dda_df = pd.DataFrame(dda_files, columns = ["Spectra_Filepath"])
    pepxml_df = pd.DataFrame(pepxml_files, columns =  ["Id_Filepath"])
    df = dda_df.reset_index().rename({"index":"Sample"}, axis = "columns")
    df["Sample"] += 1
    #df["BatchID"] = np.random.randint(1,3, size=len(df["Sample"]))    
    df["BatchID"] = 1
    df["Id_Filepath"] = pepxml_df
    df = df[headers]
    return df

mzml_dir = "dda/"
pepxml_dir = "diaumpire_spectral_lib_20210706/MSFragger_20210708_dda_only/"
#pepxml_dir = mzml_dir
spectral_lib_df = generate_input_sheet_dda(mzml_dir = mzml_dir, pepxml_dir = pepxml_dir)
spectral_lib_df.to_csv("nf_input_dda_sheet.csv", sep = "\t", index = False)
    