#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:06:50 2022

@author: ptruong
"""

import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl
sns.set_context("talk")
rcParams["text.usetex"] = True


output = "pathway_count_barplot.png"
#f, ax = plt.subplots(1, 1, figsize = (15,10))

df = pd.read_csv("pathway_count_fc_0.415_whole_bgr.tsv", sep = "\t").rename({"Unnamed: 0":"method"}, axis=1)
df["method"] = df.method.map(lambda x:x.split("_")[0])
df.rename({"method":"Method"}, axis = 1, inplace=True)

def method_renamer(x):
    if x == "reported":
        return "Reported"
    elif x == "triqler":
        return "Triqler"
    elif x == "top3":
        return "Top3"
    elif x == "msqrob2":
        return "MSqRob2"
    elif x== "msstats":
        return "MSstats"

df["Method"] = df.Method.map(lambda x:method_renamer(x))

    
df = pd.melt(df, id_vars = "Method", var_name = "group", value_name = "count")
#sns.factorplot(x='group', y='count', hue='method', data=df, kind='bar', ax= ax)
ax = sns.catplot(x='group', y='count', hue='Method', data=df, kind='bar')

#ax.set_xlabel("Cluster", fontsize=38)
#ax.set_ylabel("Number of Pathways", fontsize=38)

fp = mpl.font_manager.FontProperties(fname="/usr/share/fonts/truetype/fira/FiraSans-Regular.ttf")
plt.xlabel("Cluster", fontproperties=fp, fontsize=24,
	#labelpad=-8
)
plt.ylabel("Number of Pathways", fontproperties=fp, fontsize=24, 
#	position=(0,.3)
)

#plt.tick_params(axis='x', which='major', labelsize=40, fontproperties=fp)
#plt.tick_params(axis='y', which='major', labelsize=40, fontproperties=fp)

#ax.tick_params(axis='x', which='major', labelsize=38)#labelrotation=90)
#ax.tick_params(axis='y', which='major', labelsize=38)
#fig = plt.get_figure()
#fig.savefig(output)



