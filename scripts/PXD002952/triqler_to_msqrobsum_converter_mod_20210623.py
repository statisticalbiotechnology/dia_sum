#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 01:58:18 2021

@author: ptruong
"""

import pandas as pd

def get_mSqRobSum_input(triqler):
    def get_expr_from_triqler(triqler):
        sample_df_array = []
        samples = triqler.run.unique()
        for sample in samples:
            sample_df = triqler[triqler.run == sample]
            sample_df.sort_values("searchScore", ascending = False, inplace=True)
            sample_df = sample_df.groupby("peptide").head(1) # double check this
            run = sample_df.run.unique()[0]
            sample_df = sample_df.drop(["condition", "run", "charge", "searchScore"], axis=1).set_index(["peptide", "proteins"])
        
            sample_df[run] = sample_df["intensity"]
            sample_df.drop(["intensity"], axis=1, inplace=True)
            sample_df_array.append(sample_df)
        
        res = pd.concat(sample_df_array, axis=1)
        return res
    
    def get_fd_from_expr(res):
        fd = res.reset_index()[["peptide", "proteins"]]
        
        def map_decoy(x):
            try:
                if x.split("_")[0] == "DECOY":
                    return "TRUE"
                else:
                    return "FALSE"
            except:
                return np.nan
                
        
        map_specie = lambda x : x.split("_")[-1] 
        
        def map_specie(x):
            try:
                return x.split("_")[-1]
            except:
                return np.nan
        def map_human(x):
            if x == "HUMAN":
                return "TRUE"
            else:
                return "FALSE"
            
        def map_ecoli(x):
            if x == "ECOLI":
                return "TRUE"
            else:
                return "FALSE"
            
        def map_yeas8(x):
            if x == "YEAS8" or x == "YEAST":
                return "TRUE"
            else:
                return "FALSE"
        fd["specie"] = fd.proteins.map(map_specie)
        fd["human"] =fd.specie.map(map_human)
        fd["yeas8"] =fd.specie.map(map_yeas8)
        fd["ecoli"] =fd.specie.map(map_ecoli)
        fd.drop("specie", axis = 1, inplace = True)
        fd["decoy"]= fd.proteins.map(map_decoy)
        return fd
    
    def get_pd_from_triqler(triqler):
        pd_dat = (triqler["run"] + "_" +  triqler["condition"].astype(str)).unique()
        pd_ = pd.DataFrame(pd_dat, columns = ["sample_condition"])
        
        def sample_condition_split_ret_sample(x):
            sample = x.split("_")[0]
            return sample
        
        def sample_condition_split_ret_condition(x):
            condition = x.split("_")[1]
            return condition
        
        pd_["sample"] = pd_["sample_condition"].map(sample_condition_split_ret_sample)
        pd_["condition"] = pd_["sample_condition"].map(sample_condition_split_ret_condition)
        pd_.drop("sample_condition",axis=1, inplace = True)
        return pd_
    
    res_ = get_expr_from_triqler(triqler)
    pd_ = get_pd_from_triqler(triqler)
    fd_ = get_fd_from_expr(res_)
    expr_ = res_.reset_index().drop("proteins", axis = 1)
    return expr_, fd_, pd_

df = pd.read_csv("triqler_input_diann.csv", sep = "\t")


expr_.to_csv("expr.csv", sep = "\t")
fd_.to_csv("fd.csv", sep = "\t")

mapper = lambda x : "X" + x.split("-")[0] + "." + x.split("-")[1]
pd__ = pd_
pd__["sample"]= pd_["sample"].map(mapper)
pd_
pd__.to_csv("pd.csv", sep  = "\t")

os.("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger/msqrob_input")

/hdd_14T/data/PXD002952/20210614_dataset/diaumpire/msfragger/diann

 
et = pd.read_csv("expr.csv", sep="\t")
ft  = pd.read_csv("fd.csv", sep="\t")
pt =pd.read_csv("pd.csv", sep ="\t")

