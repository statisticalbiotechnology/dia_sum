#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:20:17 2022

@author: ptruong
"""

import pandas as pd 
import numpy as np
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *
import gseapy as gp
import io

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST



def get_uniProtKB_ID_to_KEGG_mapper(ids): # uniProtKB_ID to KEGG
    # ids = list(c1.index)
    
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db="KEGG", ids=ids
    )
    
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        # Equivalently using the stream endpoint which is more demanding
        # on the API and so is less stable:
        # results = get_id_mapping_results_stream(link)
    
    from_list = []
    to_list = []    
    for i in results["results"]:
        from_list.append(i["from"])
        to_list.append(i["to"])
    
    return dict(zip(from_list, to_list))   

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

# Some code to return a Pandas dataframe, given tabular text
def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)



def get_mouse_pathway(pathway_id = "path:mmu04210", common_name = "Apoptosis"):
    result = REST.kegg_link(target_db = "mmu", source_db = pathway_id).read()
    res = to_df(result)
    res = res.rename({0:"pathway", 1:"gene"}, axis = 1)
    res["term"] = common_name
    return res


def oxaliplatin_pathways():
    pathway_list = {"Apoptosis": "path:mmu04210",
                    "Nucleotide excision repair": "path:mmu03420",
                    "Mismatch repair": "path:mmu03430",
                    "ErbB signaling pathway": "path:mmu04012",
                    "Cell cycle": "path:mmu04110",
                    "p53 signaling pathway": "path:mmu04115"}
    
    pathway_dfs = []
    for i in pathway_list:
        pathway_dfs.append(get_mouse_pathway(pathway_list[i], common_name = i))
    res = pd.concat(pathway_dfs)
    return res    

def get_KEGG_pathway_proteins(final_df, pathway_df, pathway_term = "Apoptosis"):
    pathway_df[pathway_df.term == pathway_term]    
    res = pd.concat([final_df[final_df.KEGG_id.isin(pathway_df[pathway_df.term == pathway_term].gene)].set_index("KEGG_id"), 
               pathway_df[pathway_df.term == pathway_term].set_index("gene")], axis = 1)
    return res



os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-09-22_ctrl_vs_ST_vs_LT_study")

triqler_input = "ctrl_LT/proteins_fc_0.1"
triqler_output = triqler_input + "_pathway_ext"

triqler = parse_triqler(triqler_input)
protein_mapper = get_uniProtKB_ID_to_KEGG_mapper(triqler.protein)
triqler["KEGGProtein"] = triqler.protein.map(protein_mapper)
pathways = oxaliplatin_pathways()
for pathway in pathways.term.unique():
    triqler["pathway:" + pathway] = triqler.KEGGProtein.isin(pathways[pathways.term == pathway].gene)

triqler.to_csv(triqler_output, sep = "\t", index = False)
