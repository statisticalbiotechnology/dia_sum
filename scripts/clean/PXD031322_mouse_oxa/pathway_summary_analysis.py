#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:01:12 2022

@author: ptruong
"""


import pandas as pd
import numpy as np



reported = pd.read_csv("enrichr_reported.tsv", sep = "\t")
triqler = pd.read_csv("enrichr_triqler_fc_0.415_whole_bgr.tsv", sep = "\t")


def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

#intersection(reported[reported.group == cluster].Term, triqler[triqler.group == cluster].Term)

#reported[reported.group == cluster].append(reported[reported.group == cluster])

reported_total = pd.DataFrame()
triqler_total = pd.DataFrame()

cluster_list = []
triqler_counts = []
reported_counts = []
for cluster in triqler.group.unique():
    
    reported_count = len(reported[reported.group == cluster].Term)
    triqler_count = len(triqler[triqler.group == cluster].Term)
    
    if triqler_total.empty:
        triqler_total = triqler[triqler.group == cluster]
    else:
        triqler_total = triqler_total.append(triqler[triqler.group == cluster])
    if reported_total.empty:
        reported_total = reported[reported.group == cluster]
    else:
        reported_total = reported_total.append(reported[reported.group == cluster])
    
    print(cluster)
    print("Triqler:" + str(triqler_count))
    print("Reported:" + str(reported_count))
    triqler_counts.append(triqler_count)
    reported_counts.append(reported_count)
    cluster_list.append(cluster)
    
print("Triqler total unique pathways " + str(len(triqler_total.Term.unique())))
print("Reported total unique pathways " + str(len(reported_total.Term.unique())))

cluster_list.append("total")
triqler_counts.append(len(triqler_total.Term.unique()))
reported_counts.append(len(reported_total.Term.unique()))

pathway_count = pd.DataFrame([triqler_counts, reported_counts], columns = cluster_list, index = ["triqler_pathway_count", "reported_pathway_count"])
pathway_count.to_csv("pathway_count_fc_0.415_whole_bgr.tsv", sep = "\t")


# comparison list



pivot_reported  = reported[["Term", "Adjusted P-value", "group"]]
pivot_reported["method"] = "reported"
pivot_reported = pivot_reported.pivot(index = "Term", columns = ["group", "method"])
pivot_triqler = triqler[["Term", "Adjusted P-value", "group"]]
pivot_triqler["method"] = "triqler"
pivot_triqler = pivot_triqler.pivot(index = "Term", columns = ["group", "method"])



res = pd.concat([pivot_reported, pivot_triqler], axis = 1)
res.sort_index(axis=1, level=[0,1], ascending=[False, True], inplace = True)

res.to_csv("detailed_pathway_comparison_fc_0.415_whole_bgr.csv", sep = "\t")



