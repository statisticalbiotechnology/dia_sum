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
top3 = pd.read_csv("enrichr_top3_whole_bgr.tsv", sep = "\t")
msqrob2 = pd.read_csv("enrichr_msqrob2_whole_bgr.tsv", sep = "\t")
msstats = pd.read_csv("enrichr_msstats_whole_bgr.tsv", sep = "\t")

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

#intersection(reported[reported.group == cluster].Term, triqler[triqler.group == cluster].Term)

#reported[reported.group == cluster].append(reported[reported.group == cluster])

reported_total = pd.DataFrame()
triqler_total = pd.DataFrame()
top3_total = pd.DataFrame()
msqrob2_total = pd.DataFrame()
msstats_total = pd.DataFrame()

cluster_list = []
triqler_counts = []
reported_counts = []
top3_counts = []
msqrob2_counts = []
msstats_counts = []

for cluster in triqler.group.unique():
    
    reported_count = len(reported[reported.group == cluster].Term)
    triqler_count = len(triqler[triqler.group == cluster].Term)
    top3_count = len(top3[top3.group == cluster].Term)
    msqrob2_count = len(msqrob2[msqrob2.group == cluster].Term)
    msstats_count = len(msstats[msstats.group == cluster].Term)
    
    if triqler_total.empty:
        triqler_total = triqler[triqler.group == cluster]
    else:
        triqler_total = triqler_total.append(triqler[triqler.group == cluster])
        
    if reported_total.empty:
        reported_total = reported[reported.group == cluster]
    else:
        reported_total = reported_total.append(reported[reported.group == cluster])
        
    if top3_total.empty:
        top3_total = top3[top3.group == cluster]
    else:
        top3_total = top3_total.append(top3[top3.group == cluster])
    
    if msqrob2_total.empty:
        msqrob2_total = msqrob2[msqrob2.group == cluster]
    else:
        msqrob2_total = msqrob2_total.append(msqrob2[msqrob2.group == cluster])
        
    if msstats_total.empty:
        msstats_total = msstats[msstats.group == cluster]
    else:
        msstats_total = msstats_total.append(msstats[msstats.group == cluster])
    
    print(cluster)
    print("Triqler:" + str(triqler_count))
    print("Reported:" + str(reported_count))
    print("Top3:" + str(top3_count))
    print("Msqrob2:" + str(msqrob2_count))
    print("MSstats:" + str(msstats_count))
    triqler_counts.append(triqler_count)
    reported_counts.append(reported_count)
    top3_counts.append(top3_count)
    msqrob2_counts.append(msqrob2_count)
    msstats_counts.append(msstats_count)
    cluster_list.append(cluster)
    
print("Triqler total unique pathways " + str(len(triqler_total.Term.unique())))
print("Reported total unique pathways " + str(len(reported_total.Term.unique())))
print("Top3 total unique pathways " + str(len(top3_total.Term.unique())))
print("Msqrob2 total unique pathways " + str(len(msqrob2_total.Term.unique())))
print("MSstats total unique pathways " + str(len(msstats_total.Term.unique())))

cluster_list.append("total")
triqler_counts.append(len(triqler_total.Term.unique()))
reported_counts.append(len(reported_total.Term.unique()))
top3_counts.append(len(top3_total.Term.unique()))
msqrob2_counts.append(len(msqrob2_total.Term.unique()))
msstats_counts.append(len(msstats_total.Term.unique()))

pathway_count = pd.DataFrame([triqler_counts, reported_counts, top3_counts, msqrob2_counts, msstats_counts], columns = cluster_list, index = ["triqler_pathway_count", "reported_pathway_count", "top3_pathway_count", "msqrob2_pathway_count", "msstats_pathway_count"])
pathway_count.to_csv("pathway_count_fc_0.415_whole_bgr.tsv", sep = "\t")


# comparison list



pivot_reported  = reported[["Term", "Adjusted P-value", "group"]]
pivot_reported["method"] = "reported"
pivot_reported = pivot_reported.pivot(index = "Term", columns = ["group", "method"])
pivot_triqler = triqler[["Term", "Adjusted P-value", "group"]]
pivot_triqler["method"] = "triqler"
pivot_triqler = pivot_triqler.pivot(index = "Term", columns = ["group", "method"])

pivot_top3 = top3[["Term", "Adjusted P-value", "group"]]
pivot_top3["method"] = "top3"
pivot_top3 = pivot_top3.pivot(index = "Term", columns = ["group", "method"])
pivot_msqrob2 = msqrob2[["Term", "Adjusted P-value", "group"]]
pivot_msqrob2["method"] = "msqrob2"
pivot_msqrob2 = pivot_msqrob2.pivot(index = "Term", columns = ["group", "method"])
pivot_msstats = msstats[["Term", "Adjusted P-value", "group"]]
pivot_msstats["method"] = "msstats"
pivot_msstats = pivot_msstats.pivot(index = "Term", columns = ["group", "method"])




res = pd.concat([pivot_reported, pivot_triqler, pivot_top3, pivot_msqrob2, pivot_msstats], axis = 1)
res.sort_index(axis=1, level=[0,1], ascending=[False, True], inplace = True)

res.to_csv("detailed_pathway_comparison_fc_0.415_whole_bgr.csv", sep = "\t")



