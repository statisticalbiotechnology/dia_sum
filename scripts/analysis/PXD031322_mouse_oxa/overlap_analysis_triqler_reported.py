#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:04:55 2022

@author: ptruong
"""

import pandas as pd 
import numpy as np
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *


def get_mapped_proteins(ids):
    # ids = list(c1.index)
    
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=ids
    )
    
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        # Equivalently using the stream endpoint which is more demanding
        # on the API and so is less stable:
        # results = get_id_mapping_results_stream(link)
    
    
    primaryAccession_list = []
    secondaryAccessions_list = []
    uniProtKB_ID_list = []
    geneName_list = []
    geneName_synonyms_list = []
    
    for i in range(len(results["results"])):    
        #print(i)
        primaryAccession = results["results"][i]["to"]["primaryAccession"]
        try:
            secondaryAccessions = ";".join(results["results"][i]["to"]["secondaryAccessions"])
        except:
            secondaryAccessions = ""
        uniProtKB_ID = results["results"][i]["to"]["uniProtkbId"]
        try:
            geneName = results["results"][i]["to"]["genes"][0]["geneName"]["value"]
        except:
            geneName = ""
        synonyms_list = []
        try:
            for s in results["results"][i]["to"]["genes"][0]["synonyms"]:
                synonyms_list.append(s["value"])
        except:
            pass
        geneName_synonyms = ";".join(synonyms_list)
        
        primaryAccession_list.append(primaryAccession)
        secondaryAccessions_list.append(secondaryAccessions)
        uniProtKB_ID_list.append(uniProtKB_ID)
        geneName_list.append(geneName)
        geneName_synonyms_list.append(geneName_synonyms)
    
    res = pd.DataFrame([primaryAccession_list, secondaryAccessions_list,
                  uniProtKB_ID_list, geneName_list, geneName_synonyms_list],
                 index = ["primaryAccession", "secondaryAccession", "uniProtKB_ID", "geneName", "geneNameSynonyms"]).T
    return res    


def intersection(a,b):
    return list(set(a)&(set(b)))


def intersection_count(triqler, reported):
    groups = ["C1", "C2", "C3", "C4", "C5", "C6"]
    
    triqler_protein_count_list = []
    reported_protein_count_list = []
    intersection_protein_count_list = []
    diff_triqler_count_list = []
    diff_reported_count_list = []
    for group in groups:
        triqler_group = triqler[triqler["Triqler"] == group]
        reported_group = reported[reported["Reported"] == group]
        intersecting_proteins = intersection(triqler_group.protein,reported_group.protein)
        
        count_triqler_proteins = len(triqler_group)
        count_reported_proteins = len(reported_group)
        count_intersecting_proteins = len(intersecting_proteins)
        diff_triqler = count_triqler_proteins - count_intersecting_proteins
        diff_reported = count_reported_proteins - count_intersecting_proteins
    
        triqler_protein_count_list.append(count_triqler_proteins)
        reported_protein_count_list.append(count_reported_proteins)
        intersection_protein_count_list.append(count_intersecting_proteins)
        diff_triqler_count_list.append(diff_triqler)
        diff_reported_count_list.append(diff_reported)
    
    res = pd.DataFrame([triqler_protein_count_list, reported_protein_count_list, intersection_protein_count_list,
                  diff_triqler_count_list, diff_reported_count_list], 
                 columns = groups,
                 index = ["triqler_count", "reported_count", "intersection_count",
                          "triqler-intersection", "reported-intersection"]).T    
    return res    



reported = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
reported = reported[~reported.Condition.isin(["Others"])]
reported = reported[["Protein.Ids", "Condition"]].rename({"Protein.Ids":"protein", "Condition":"Reported"},axis=1)
mapped_proteins = get_mapped_proteins(list(reported.protein))
mapper = mapped_proteins[["primaryAccession", "uniProtKB_ID"]].set_index("primaryAccession").to_dict()["uniProtKB_ID"]
reported["protein"] = reported.protein.map(mapper)

triqler = pd.read_csv("triqler_protein_fc_0.415.tsv", sep = "\t").rename({"0":"Triqler", "Unnamed: 0":"protein"}, axis = 1)
#mapped_proteins = get_mapped_proteins(list(triqler.protein))
#mapper = mapped_proteins[["uniProtKB_ID", "primaryAccession"]].set_index("uniProtKB_ID").to_dict()["primaryAccession"]
#triqler["mapped_protein"] = triqler.protein.map(mapper)

res = intersection_count(triqler = triqler, reported = reported)
res.to_csv("triqler_vs_reported_intersection.tsv", sep = "\t")





