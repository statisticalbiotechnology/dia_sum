#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:14:14 2022

@author: ptruong
"""

import pandas as pd
import numpy as np 

df = pd.read_csv("PBMC_5000_DIA-NN_GPF_based_report.tsv", sep = "\t")

df.columns
df["Sample"] = df.Run.map(lambda x:x.split("_")[0])
df["Replicate"] = df.Run.map(lambda x:x.split("_")[1])
df.drop_duplicates(inplace = True)

df[df["Protein.Group"] == "Q9Y4H2"].to_csv("checkthis.tsv", sep = "\t", index = False)
df["Precursor.Quantity/Charge"] = df["Precursor.Quantity"] / df["Precursor.Charge"]


df_subset = df[["Protein.Names", "Stripped.Sequence", "Protein.Group", "Run", "Precursor.Quantity/Charge", "Precursor.Quantity", "Precursor.Charge"]]
df_subset.drop_duplicates(inplace=True)
df_s = df_subset[df_subset["Protein.Group"] == "Q9Y4H2"]
df_s.drop_duplicates()

df_s.groupby(["Protein.Names", "Run"])["Precursor.Quantity/Charge"].max()


df_s.pivot(index = ["Protein.Names", "Protein.Group", "Stripped.Sequence", "Precursor.Charge"], columns = "Run")
df_s

pivot = df_subset.pivot(index = ["Protein.Group", "Precursor.Charge"], columns = "Run", values = "Precursor.Quantity")


report = pd.read_csv("report-first-pass.tsv", sep = "\t")
report[report["Protein.Group"] == "P55011"].to_csv("P55011.tsv", sep = "\t", index = False)

pg_matrix = pd.read_csv("report-first-pass.pg_matrix.tsv", sep = "\t")
pg_matrix[pg_matrix["Protein.Group"] == "P55011"].T

#df[df["Protein.Group"] == "P55011"]


# Work from here to see number of diff. exp genes.
import qnorm

pg_matrix.set_index("Protein.Group", inplace = True)
pg_matrix = pg_matrix.iloc[:,pg_matrix.columns.str.contains("mzML")]
#qnorm.quantile_normalize(pg_matrix, axis = 1)
pg_matrix = pg_matrix[pg_matrix.isna().sum(axis = 1) < 4] # filter away proteins with more than 50% missing values
qnorm_pg_matrix = qnorm.quantile_normalize(pg_matrix, axis = 1)
log2_qnorm_pg_matrix = np.log2(qnorm_pg_matrix)
#qnorm.quntile_normalize seem to produce some  nans
log2_qnorm_pg_matrix = log2_qnorm_pg_matrix[log2_qnorm_pg_matrix.isna().sum(axis = 1) < 4]
log2_qnorm_pg_matrix.to_csv("log2_qnorm_pg_matrix.tsv", sep = "\t", index = True)

#pg_matrix[pg_matrix.index == "P62805"]
#qnorm_pg_matrix[qnorm_pg_matrix.index == "P62805"]
#log2_qnorm_pg_matrix[log2_qnorm_pg_matrix.index == "P62805"]

# R for imputeLCMD
imputed_log2_qnorm_pg_matrix = pd.read_csv("imputed_log2_qnorm_pg_matrix.tsv", sep = "\t")
col = imputed_log2_qnorm_pg_matrix.columns
multiCol = pd.MultiIndex.from_tuples([(col.split("_")[-2] + "_" + col.split("_")[-1].split(".")[0], 
   col.split("_")[-2], 
   col.split("_")[-1].split(".")[0]) for col in imputed_log2_qnorm_pg_matrix], 
 names = ["sample_replicate", "sample", "replicate"])
imputed_log2_qnorm_pg_matrix.columns = multiCol

df = imputed_log2_qnorm_pg_matrix
ctrl = df.iloc[:,df.columns.get_level_values("sample") == "Control"].mean(axis = 1)
case = df.iloc[:,df.columns.get_level_values("sample") == "LPS"].mean(axis = 1)


res = pd.DataFrame(ctrl - case)
res[abs(res) > 2].dropna() # 573 differentially expressed proteins.







