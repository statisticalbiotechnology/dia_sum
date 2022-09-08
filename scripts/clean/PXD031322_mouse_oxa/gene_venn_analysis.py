#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 10:53:27 2022

@author: ptruong
"""



import pandas as pd 
import numpy as np
import os
os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322")


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

def get_clusters(ctrl_lt_file = "proteins.1vs2.tsv",
                 ctrl_st_file = "proteins.1vs3.tsv",
                 lt_st_file = "proteins.2vs3.tsv"):
    cols = ['q_value', 'log2_fold_change', "upregulated"]
    
    ctrl_lt = parse_triqler(ctrl_lt_file).set_index("protein")
    ctrl_st = parse_triqler(ctrl_st_file).set_index("protein")
    lt_st = parse_triqler(lt_st_file).set_index("protein")
    
    ctrl_lt["upregulated"] = ctrl_lt.log2_fold_change > 0
    ctrl_st["upregulated"] = ctrl_st.log2_fold_change > 0 
    lt_st["upregulated"] = lt_st.log2_fold_change > 0
    
    
    ctrl_lt = ctrl_lt[cols].rename({"q_value":"ctrl-lt:q_value", "upregulated":"ctrl-lt:upregulated", "log2_fold_change":"ctrl-lt:log2_fold_change"}, axis = 1)
    ctrl_st = ctrl_st[cols].rename({"q_value":"ctrl-st:q_value", "upregulated":"ctrl-st:upregulated", "log2_fold_change":"ctrl-st:log2_fold_change"}, axis = 1)
    lt_st = lt_st[cols].rename({"q_value":"lt-st:q_value", "upregulated":"lt-st:upregulated", "log2_fold_change":"lt-st:log2_fold_change"}, axis = 1)
    
    df = pd.concat([ctrl_lt, ctrl_st, lt_st], axis = 1)
    
    # Clusters 
    C1 = ((df["ctrl-lt:upregulated"] == True) & (df["ctrl-st:upregulated"] == True) & (df["lt-st:upregulated"] == False))
    C2 = ((df["ctrl-lt:upregulated"] == False) & (df["ctrl-st:upregulated"] == False) & (df["lt-st:upregulated"] == True))
    C3 = (df["lt-st:upregulated"] == False)      
    C4 = (df["lt-st:upregulated"] == True)
    C5 = (df["ctrl-st:upregulated"] == True)
    C6 = (df["ctrl-st:upregulated"] == False)
    return C1, C2, C3, C4, C5, C6

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


df = pd.read_excel("reported_6_subclusters.xlsx", header = 1)

reported_clusters = {}
for cluster in df[df.Condition != "Others"].Condition.unique():
    
    df_cluster = df[df.Condition == cluster]
    
    print(cluster)
    print("ST-CTRL p-value max" + str(df_cluster["ST-Ctrl_p adj"].max()))
    print("LT-CTRL p-value max" + str(df_cluster["LT-Ctrl_p adj"].max()))
    print("LT-ST p-value max" + str(df_cluster["LT-ST_p adj"].max()))
    print("")
    
    reported_clusters[cluster] = df[df.Condition == cluster]
    
        
os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-08-11_run/2022-09-09_triqler_and_enrichr_results/fc_0.415")

ctrl_lt = parse_triqler("proteins.1vs2.tsv").set_index("protein")
ctrl_st = parse_triqler("proteins.1vs3.tsv").set_index("protein")
lt_st = parse_triqler("proteins.2vs3.tsv").set_index("protein")
C1, C2, C3, C4, C5, C6 = get_clusters(ctrl_lt_file = "proteins.1vs2.tsv",
                                      ctrl_st_file = "proteins.1vs3.tsv",
                                      lt_st_file = "proteins.2vs3.tsv")

fdr_threshold = 0.05

bgr_proteins = pd.concat([ctrl_lt.reset_index(),ctrl_st.reset_index(), lt_st.reset_index()]).protein.unique()

# all proteins

# Check ctrl-lt C1,C2
c1 = ctrl_lt[C1 & (ctrl_lt.q_value < fdr_threshold)]
c2 = ctrl_lt[C2 & (ctrl_lt.q_value < fdr_threshold)]

# Check lt-st C3, c4
c3 = lt_st[C3 & (lt_st.q_value < fdr_threshold)]
c4 = lt_st[C4 & (lt_st.q_value < fdr_threshold)]

# Check ctrl-st C5, C6
c5 = ctrl_lt[C5 & (ctrl_lt.q_value < fdr_threshold)]
c6 = ctrl_lt[C6 & (ctrl_lt.q_value < fdr_threshold)]


c1_mapped = get_mapped_proteins(ids = list(c1.index))
c2_mapped = get_mapped_proteins(ids = list(c2.index))
c3_mapped = get_mapped_proteins(ids = list(c3.index))
c4_mapped = get_mapped_proteins(ids = list(c4.index))
c5_mapped = get_mapped_proteins(ids = list(c5.index))
c6_mapped = get_mapped_proteins(ids = list(c6.index))


triqler_clusters = dict(zip([f"C{i+1}" for i in range(6)], [c1_mapped, c2_mapped, c3_mapped, c4_mapped, c5_mapped, c6_mapped]))


cluster = "C1"

triqler_cluster = triqler_clusters[cluster] #We need to check synonyms as well
# If gene matches main gene or one of the synonyms, its an overlap...
# So we need to for-loop

triqler_clusters[cluster].columns
triqler_cluster["genes"] = triqler_cluster.geneName +";" +  triqler_cluster.geneNameSynonyms

overlap = 0
for reported_gene in reported_clusters[cluster]["Gene.Symbol"]:
    for triqler_gene_and_geneSynonyms in triqler_cluster.genes:
        triqler_match = False
        for triqler_gene in triqler_gene_and_geneSynonyms.split(";"):
            if reported_gene == triqler_gene:
                overlap += 1
                triqler_match = True
                break
            print(triqler_gene + ":" + reported_gene)
        if triqler_match == True:
            triqler_match = False
            break
    
            
    











