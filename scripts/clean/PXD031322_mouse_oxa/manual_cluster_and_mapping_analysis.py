#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:54:35 2022

@author: ptruong

https://www.uniprot.org/help/id_mapping
https://biit.cs.ut.ee/gprofiler/page/apis


ToDo:
    
    Continue on gprofiler background
    Plot gprofiler KEGG pathways as in paper
"""
import pandas as pd 
import numpy as np
from uniprot_idmapper import *


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

ctrl_lt = parse_triqler("proteins.1vs2.tsv").set_index("protein")
ctrl_st = parse_triqler("proteins.1vs3.tsv").set_index("protein")
lt_st = parse_triqler("proteins.2vs3.tsv").set_index("protein")


fdr_threshold = 0.05
# Check ctrl-lt C1,C2
c1 = ctrl_lt[C1 & (ctrl_lt.q_value < fdr_threshold)]
c2 = ctrl_lt[C2 & (ctrl_lt.q_value < fdr_threshold)]

# Check lt-st C3, c4
c3 = lt_st[C1 & (lt_st.q_value < fdr_threshold)]
c4 = lt_st[C2 & (lt_st.q_value < fdr_threshold)]

# Check ctrl-st C5, C6
c5 = ctrl_lt[C1 & (ctrl_lt.q_value < fdr_threshold)]
c6 = ctrl_lt[C2 & (ctrl_lt.q_value < fdr_threshold)]


c1.index

c1_mapped = get_mapped_proteins(ids = list(c1.index))
c2_mapped = get_mapped_proteins(ids = list(c2.index))
c3_mapped = get_mapped_proteins(ids = list(c3.index))
c4_mapped = get_mapped_proteins(ids = list(c4.index))
c5_mapped = get_mapped_proteins(ids = list(c5.index))
c6_mapped = get_mapped_proteins(ids = list(c6.index))


c1_mapped
c2_mapped
c3_mapped
c4_mapped
c5_mapped
c6_mapped

query_genes = list(c1_mapped.geneName)
user_threshold = 0.05
import requests
r = requests.post(
    url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
    json={
    'organism':'mmusculus',
    'query':query_genes,
    'sources' :['KEGG'], #only look into Gene Ontology terms.
    'user_threshold':user_threshold, #reduce the significance threshold,
    'significance_threshold_method':'bonferroni', #use bonferroni correction instrad of the default 'g_SCS'.
    'no_evidences':True, #skip lookup for evidence codes. Speeds up queries, if there is no interest in evidence codes.
    'no_iea':True, #Ignore electonically annotated GO annotations

    'domain_scope':'custom',#use the genes in the probe as the statistical background.
    'background':'AFFY_HG_U133A'
    },
    headers={
    'User-Agent':'FullPythonRequest'
    }
)

obj = r.json()['result']

type(obj[0])

