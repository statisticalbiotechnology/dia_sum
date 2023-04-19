#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 09:36:40 2021

@author: ptruong
"""

import os
import pandas as pd
import numpy as np

os.chdir("/home/ptruong/git/dia_sum/scripts/PXD002952")
from triqler_output_to_df import parse_triqler

os.chdir("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix/triqler_human_fc05_full_tsv")



df = parse_triqler("proteins.tsv")

df = df[df.q_value <  0.05]


df = parse_triqler("proteins.tsv")
#df = df[df.protein_id_posterior_error_prob < 0.01]
specie_mapper = lambda x : x.split("_")[-1]
df["specie"] = df.protein.map(specie_mapper) 


df[df.specie == "ECOLI"].tail(50).protein.to_csv("protein_list_ECOLI_n50_tail.txt", sep = "\t", index = False, header = False)
df[df.specie == "HUMAN"].tail(50).protein.to_csv("protein_list_HUMAN_n50_tail.txt", sep = "\t", index = False, header = False)
df[df.specie == "YEAS8"].tail(50).protein.to_csv("protein_list_YEAS8_n50_tail.txt", sep = "\t", index = False, header = False)

df[df.specie == "ECOLI"].sample(50).protein.to_csv("protein_list_ECOLI_n50_random_nonSig.txt", sep = "\t", index = False, header = False)
df[df.specie == "HUMAN"].tail(50).protein.to_csv("protein_list_HUMAN_n50_tail.txt", sep = "\t", index = False, header = False)
df[df.specie == "YEAS8"].sample(50).protein.to_csv("protein_list_YEAS8_n50_random_nonSig.txt", sep = "\t", index = False, header = False)


df[df.specie == "HUMAN"].to_csv("human_proteins.tsv", sep ="\t", index = False)



df

def pq_plot(df, title = "title", label = "TEST"):
    """
    Input - untresholded triqler proteins.tsv 
    """
    q_vals = []
    n_de = []
    for q in np.arange(0,1,0.01):
        q_vals.append(q)
        n_de.append(len(df[df.q_value <  q]))
    
    df_pq = pd.DataFrame(np.array([q_vals, n_de]).T, columns = ["q_value", "n_de"])
    
    import matplotlib.pyplot as plt
    
    plt.plot(df_pq.q_value, df_pq.n_de, label = label)
    plt.title(title)
    #plt.set_label(label)
#    line, = ax.plot(x, y, 'b.-', ...)
#line.set_label('line 1')

pq_plot(df, "all, --fold_change_eval = 0.5", "all")
pq_plot(df[df.specie == "YEAS8"], "YEAS8, --fold_change_eval = 0.5", "YEAS8")
pq_plot(df[df.specie == "ECOLI"], "ECOLI, --fold_change_eval = 0.5", "ECOLI")
pq_plot(df[df.specie == "HUMAN"], "HUMAN, --fold_change_eval = 0.5", "HUMAN")




