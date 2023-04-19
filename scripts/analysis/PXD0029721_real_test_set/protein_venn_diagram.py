#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:01:13 2022

@author: ptruong
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def parse_triqler(triqler_output_file):
    """
    Parses triqler output format to pandas dataframe.
    """
    f = open(triqler_output_file, "r")
    lines = f.readlines()
    line = lines.pop(0)
    cols = line.split("\n")[0].split("\t")[:]
    n_cols = len(cols)
    
    data_array = []
    for line in lines:
        line = line.split("\n")[0].split("\t")
        vals = line[:n_cols-1]
        peptides = ";".join(line[n_cols-1:])
        data = vals + [peptides]
        data_array.append(data)
    df = pd.DataFrame(data_array, columns = cols)

    df = pd.concat([df[["protein", "peptides"]], df.drop(["protein", "peptides"], axis = 1).astype(float)], axis = 1)
    
    return df


reported = pd.read_excel("1-s2.0-S0731708522002163-mmc6.xlsx")
triqler = parse_triqler("proteins.tsv")
triqler = triqler[~triqler.protein.str.contains("DECOY")]
#triqler[triqler.protein_id_posterior_error_prob < 0.05]
triqler = triqler[triqler.protein_id_posterior_error_prob < 0.01]

reported_proteins = reported.SYMBOL
triqler_proteins = triqler.protein.map(lambda x:x.split("_")[0])

intersection = list(set(reported_proteins) & set(triqler_proteins))

def plot_venn(a, b, intersect, a_name = "reported", b_name = "Triqler"):
    venn2(subsets = (len(a)-len(intersect), len(b)-len(intersect),
                     len(intersect)), set_labels = (a_name, b_name))
    plt.show()

# plot venn diagram of identified proteins
plot_venn(reported_proteins, triqler_proteins, intersection, a_name = "reported", b_name = "Triqler")


reported_significant = reported[reported.Significance.isin(["down", "up"])].SYMBOL
triqler_significant = triqler[triqler.q_value < 0.05].protein.map(lambda x:x.split("_")[0])
intersection_significant = list(set(reported_significant) & set(triqler_significant))
plot_venn(reported_significant, triqler_significant, intersection_significant, a_name = "reported", b_name = "Triqler")

reported_significant.to_csv("reported_significant.txt", index = False)
triqler_significant.to_csv("triqler_significant.txt", index = False)
pd.Series(intersection_significant).to_csv("intersection_significant.txt", index = False)

