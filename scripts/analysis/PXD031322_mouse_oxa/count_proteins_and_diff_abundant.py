#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 10:47:39 2022

@author: ptruong
"""


import pandas as pd
import numpy as np
import seaborn as sns
import os


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

def count_proteins(fdr_threshold):
    fdr_threshold = fdr_threshold
    n_proteins_list = []
    input_fdr_threshold_list = []
    for file in os.listdir():
            print("parsing file:" + file)
            res = parse_triqler(file)
            n_proteins = (res.protein_id_posterior_error_prob < fdr_threshold).sum()
            input_fdr_threshold = float(file.split("_")[-1][:-4])
            n_proteins_list.append(n_proteins)
            input_fdr_threshold_list.append(input_fdr_threshold)
    output = pd.DataFrame([n_proteins_list, input_fdr_threshold_list], index = ["n_proteins", "triqler_input_fdr_threshold_level"]).T
    output = output.sort_values(by = "triqler_input_fdr_threshold_level").reset_index().drop("index", axis = 1)
    return output

def count_diff_abundant(fdr_threshold):
    fdr_threshold = fdr_threshold
    n_diff_abundant_list = []
    input_fdr_threshold_list = []
    for file in os.listdir():
        print("parsing file:" + file)
        res = parse_triqler(file)
        n_diff_abundant = (res.q_value < fdr_threshold).sum()
        input_fdr_threshold = float(file.split("_")[-1][:-4])
        n_diff_abundant_list.append(n_diff_abundant)
        input_fdr_threshold_list.append(input_fdr_threshold)
    output = pd.DataFrame([n_diff_abundant_list, input_fdr_threshold_list], index = ["n_differentially_abundant", "triqler_input_fdr_threshold_level"]).T
    output = output.sort_values(by = "triqler_input_fdr_threshold_level").reset_index().drop("index", axis = 1)
    return output


fc_eval = 0.4
proteins = count_proteins(0.05)
title = f"Significant proteins at triqler fc_eval={str(fc_eval)} and protein_id_PEP < 0.05"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_proteins", data = proteins).set(title = title)

proteins = count_proteins(0.01)
title = f"Significant proteins at triqler fc_eval={str(fc_eval)} and protein_id_PEP < 0.01"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_proteins", data = proteins).set(title = title)

proteins = count_proteins(0.001)
title = f"Significant proteins at triqler fc_eval={str(fc_eval)} and protein_id_PEP < 0.01"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_proteins", data = proteins).set(title = title)


diff_abundant = count_diff_abundant(0.05)
title = f"Differentially abundant proteins at triqler fc_eval={str(fc_eval)} and q_value < 0.05"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_differentially_abundant", data = diff_abundant).set(title = title)

diff_abundant = count_diff_abundant(0.01)
title = f"Differentially abundant proteins at triqler fc_eval={str(fc_eval)} and q_value < 0.01"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_differentially_abundant", data = diff_abundant).set(title = title)

diff_abundant = count_diff_abundant(0.001)
title = f"Differentially abundant proteins at triqler fc_eval={str(fc_eval)} and q_value < 0.001"
sns.lineplot(x = "triqler_input_fdr_threshold_level", y = "n_differentially_abundant", data = diff_abundant).set(title = title)













