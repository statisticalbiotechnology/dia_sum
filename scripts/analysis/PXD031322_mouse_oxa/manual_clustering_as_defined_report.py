#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 11:28:15 2022

@author: ptruong
"""


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

# pd.DataFrame(c1.index).to_csv("c1_fdr_0.05_fc_eval_0.6.tsv", header = False, index = False)
c1.index


pd.DataFrame(c1.index).to_csv("test_rows.csv", sep = ",", index = False, header = False)
pd.DataFrame(c1.index).T.to_csv("test.csv", sep = ",", index = False, header = False)



# Get gene names
###########


import requests, sys
from urllib import request
from bs4 import BeautifulSoup
import json



accession = "P21802"
def get_gene_name_from_protein_accession(accession = "P21802"):
#    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession={accession}"
    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&id={accession}"
    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&gene=TRFL"
    r = requests.get(requestURL, headers={"Accept":"application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    responseBody = r.text
    
    
    obj = json.loads(responseBody)
    return obj[0]["gene"][0]["name"]["value"]


--form 'from=UniProtKB_AC-ID' \
     --form 'to=UniProtKB' \
         
         
         
for key in obj[0]:
    print(key)

obj[0]["id"]

get_gene_name_from_protein_accession(accession = "P21802")

c1.reset_index().protein.map(lambda x:x.split("_")[0])[0]




get_gene_name_from_protein_accession(accession = "TRFL")

# https://www.uniprot.org/tool-dashboard, use this one for mappinh


genes = pd.read_csv("uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.08.30-06.58.27.06.tsv", sep = "\t")
genes["first_gene_name"] = genes["Gene Names"].map(lambda x:x.split(" ")[0])
genes["first_gene_name"].to_csv("genes_c1_fdr_0.05_fc_eval_0.6.tsv", index = False)























