#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:40:11 2023

@author: ptruong
"""

import pandas as pd
from parsers.parse_triqler import *
import os

os.chdir("/home/ptruong/git/dia_sum/scripts/clean/results/PS")
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/results/ID")
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa/data")

top_n = 10 
msstats = pd.read_csv("msstats_adj_results.csv", sep = ",")
#msstats = pd.read_csv("msstats_results.csv", sep = ",") #LFQBench
msstats = msstats.sort_values(by = "adj.pvalue").reset_index().drop("index", axis = 1)

msstats_ranking = pd.DataFrame(
    pd.concat([
        pd.DataFrame(msqrob2[msqrob2.protein.isin(msstats.head(top_n).Protein)].index), 
        pd.DataFrame(top3[top3.ProteinName.isin((msstats.head(top_n).Protein))].index),
        pd.DataFrame(triqler[triqler.protein.isin((msstats.head(top_n).Protein))].index),
        pd.DataFrame(msstats.head(top_n).Protein)], axis = 1).values, columns = ["msqrob2","top3", "triqler", "msstats_top_protein"]
             ).set_index("msstats_top_protein")

msqrob2 = pd.read_csv("msqrob2_results.csv", sep = ",").rename({"Unnamed: 0": "protein"}, axis = 1)
msqrob2 = msqrob2.sort_values(by = "adjPval").reset_index().drop("index", axis = 1)
msqrob2_ranking = pd.DataFrame(
    pd.concat([
        pd.DataFrame(msstats[msstats.Protein.isin(msqrob2.head(top_n).protein)].index), 
        pd.DataFrame(top3[top3.ProteinName.isin(msqrob2.head(top_n).protein)].index),
        pd.DataFrame(triqler[triqler.protein.isin(msqrob2.head(top_n).protein)].index),
        pd.DataFrame(msqrob2.head(top_n).protein)
        
        
        ], axis = 1).values, columns = ["msstats","top3", "triqler", "msqrob2_top_protein"]
             ).set_index("msqrob2_top_protein")


msqrob2[msqrob2.protein.isin(top3.head(top_n).ProteinName)]
msqrob2[msqrob2.protein.isin(triqler.head(top_n).protein)]

top3 = pd.read_csv("top3_results.csv", sep = "\t")
top3 = top3.sort_values(by = "q").reset_index().drop("index", axis = 1)
top3_ranking = pd.DataFrame(
    pd.concat([
        pd.DataFrame(msstats[msstats.Protein.isin(top3.head(top_n).ProteinName)].index), 
        pd.DataFrame(msqrob2[msqrob2.protein.isin(top3.head(top_n).ProteinName)].index),
        pd.DataFrame(triqler[triqler.protein.isin(top3.head(top_n).ProteinName)].index),
        pd.DataFrame(top3.head(top_n).ProteinName)], axis = 1).values, columns = ["msstats","msqrob2", "triqler", "top3_top_protein"]
             ).set_index("top3_top_protein")



triqler = parse_triqler("triqler_results.tsv")
#triqler = parse_triqler("triqler_results.csv") # LFQBench
triqler_ranking = pd.DataFrame(
    pd.concat([
        pd.DataFrame(msstats[msstats.Protein.isin(triqler.head(top_n).protein)].index), 
        pd.DataFrame(msqrob2[msqrob2.protein.isin(triqler.head(top_n).protein)].index),
        pd.DataFrame(top3[top3.ProteinName.isin(triqler.head(top_n).protein)].index),
        pd.DataFrame(triqler.head(top_n).protein)], axis = 1).values, columns = ["msstats","msqrob2", "top3", "triqler_top_protein"]
             ).set_index("triqler_top_protein")

msstats_ranking.to_csv("msstats_ranking.tsv", sep = "\t")
msqrob2_ranking.to_csv("msqrob2_ranking.tsv", sep = "\t")
top3_ranking.to_csv("top3_ranking.tsv", sep = "\t")
triqler_ranking.to_csv("triqler_ranking.tsv", sep = "\t")


