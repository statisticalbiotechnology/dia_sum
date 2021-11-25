#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 23:45:17 2021

@author: ptruong
"""

import pandas as pd
import os 
import matplotlib.pyplot as plt


os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")


df1 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_002-Pedro_-_Sample_1_-_SW32_-_Repl1.mzML_with_dscore.csv", sep = "\t")
df2 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_003-Pedro_-_Sample_2_-_SW32_-_Repl1.mzML_with_dscore.csv", sep = "\t")
df3 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_004-Pedro_-_Sample_1_-_SW32_-_Repl2.mzML_with_dscore.csv", sep = "\t")
df4 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_005-Pedro_-_Sample_2_-_SW32_-_Repl2.mzML_with_dscore.csv", sep = "\t")
df5 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_006-Pedro_-_Sample_1_-_SW32_-_Repl3.mzML_with_dscore.csv", sep = "\t")
df6 = pd.read_csv("osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_007-Pedro_-_Sample_2_-_SW32_-_Repl3.mzML_with_dscore.csv", sep = "\t")


def return_sorted_m_score(df):
    df_sorted = df.sort_values(by = "m_score").m_score
    df_m = df_sorted.reset_index().drop("index", axis = 1)
    return df_m


df1_m = return_sorted_m_score(df1)
df2_m = return_sorted_m_score(df2)
df3_m = return_sorted_m_score(df3)
df4_m = return_sorted_m_score(df4)
df5_m = return_sorted_m_score(df5)
df6_m = return_sorted_m_score(df6)

df1_m["002_sample1"] = df1_m["m_score"]
df2_m["003_sample2"] = df2_m["m_score"]
df3_m["004_sample1"] = df3_m["m_score"]
df4_m["005_sample2"] = df4_m["m_score"]
df5_m["006_sample1"] = df5_m["m_score"]
df6_m["007_sample2"] = df6_m["m_score"]

df_concat = pd.concat([df1_m["002_sample1"], df2_m["003_sample2"], df3_m["004_sample1"], df4_m["005_sample2"], df5_m["006_sample1"], df6_m["007_sample2"]], axis = 1)


ax = df_concat.plot()
ax.set_xlabel("ordered PSM by m_Score")
ax.set_ylabel("m_score")










