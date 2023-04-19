#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:59:02 2022

@author: ptruong
"""

import pandas as pd

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


df_orig = parse_triqler("test_orignal_popt.csv")
df_new = parse_triqler("test (copy).csv")


intersecting_proteins = list(set(df_orig[df_orig.q_value < 0.05].protein) & set(df_new[df_new.q_value < 0.05].protein))

df_orig = df_orig[df_orig.q_value < 0.05]
prot = df_orig[~df_orig.protein.isin(intersecting_proteins)] # are the proteins in top3, msstats, msqrob

msstats = pd.read_csv("msstats_results.csv", sep = ",")
msstats[msstats.Protein.isin(prot)]

msqrob = pd.read_csv("msqrob2_results.csv", sep = ",").rename({"Unnamed: 0": "protein"}, axis = 1)
msqrob[msqrob.protein.isin(prot)]

top3 = pd.read_csv("top3_results.csv", sep = "\t")
top3[top3.ProteinName.isin(prot)]

# check diff in log2FC
df_orig = parse_triqler("test_orignal_popt.csv")
df_new = parse_triqler("test (copy).csv")

df_orig.set_index("protein", inplace = True)
df_new.set_index("protein", inplace = True)
df_new["new_log2fc"] = df_new.log2_fold_change
df_orig["orig_log2fc"] = df_orig.log2_fold_change

df_new["new_q_value"] = df_new.q_value
df_orig["orig_q_value"] = df_orig.q_value

diff = pd.concat([df_new.new_log2fc, df_orig.orig_log2fc, df_new.new_q_value, df_orig.orig_q_value], axis = 1)
diff["diff"]  = diff.new_log2fc - diff.orig_log2fc
diff["diff_q"]  = diff.new_q_value - diff.orig_q_value

diff["diff"].max()
diff["diff"].min()
diff["diff_q"].max()
diff["diff_q"].min()

diff["specie"] = diff.index.map(lambda x:x.split("_")[-1])

diff.to_csv("diff.tsv", sep = "\t")
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "HUMAN")]
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "HUMAN")].describe()["diff"]
#diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "HUMAN")].count()
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "ECOLI")]
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "ECOLI")].describe()["diff"]
#diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "ECOLI")].count()
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "YEAST")]
diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "YEAST")].describe()["diff"]
#diff[(abs(diff["diff"]) > 0.05) & (diff.specie == "YEAST")].count()




import seaborn as sns
sns.displot(diff, x = "diff")


diff.max()
diff.min()

diff
df_orig.columns

