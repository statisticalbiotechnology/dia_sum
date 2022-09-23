#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:50:55 2022

@author: ptruong
"""


import pandas as pd 
import numpy as np
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *
import gseapy as gp


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


def intersection(a,b):
    return list(set(a)&(set(b)))


def S3_statistics(reported_S3):
    print("Printing statistics on S3 reported supplementary file.")
    print("Number of significant according to reported: " +str(len(reported_S3[~reported_S3["Condition"].isin(["Not_sig"])])))
    print("Adj.p_value max: " + str(reported_S3[~reported_S3["Condition"].isin(["Not_sig"])]["Adjusted_P_value"].max()))
    print("Adj.p_value min: " + str(reported_S3[~reported_S3["Condition"].isin(["Not_sig"])]["Adjusted_P_value"].min()))
    print("log2 (FC) max: " + str(reported_S3[~reported_S3["Condition"].isin(["Not_sig"])]['log2（FC）'].max()))
    print("log2 (FC) min: " + str(reported_S3[~reported_S3["Condition"].isin(["Not_sig"])]['log2（FC）'].min()))
    print("log2 (FC) min on Up condition: " + str(reported_S3[reported_S3["Condition"].isin(["Up"])]['log2（FC）'].min())) # larger than 0.6
    print("log2 (FC) max on Down condition: " + str(reported_S3[reported_S3["Condition"].isin(["Down"])]['log2（FC）'].max())) # smaller than -0.6

def S4_statistics(reported_S4):
    print("Printing statistics on S4 reported supplmentary file.")
    print("Adj.p_value max: " + str(reported_S4[reported_S4["Condition"].isin(["C1", "C2"])]["LT-Ctrl_p adj"].max()))
    print("Adj.p_value min: " + str(reported_S4[reported_S4["Condition"].isin(["C1", "C2"])]["LT-Ctrl_p adj"].min()))
    print("C1 min log2(FC) : " + str(reported_S4[reported_S4["Condition"].isin(["C1"])]["LT-Ctrl_diff"].min()))
    print("C2 max log2(FC) : " + str(reported_S4[reported_S4["Condition"].isin(["C2"])]["LT-Ctrl_diff"].max()))

def get_mapper(reported, protein_col):
    mapped_proteins = get_mapped_proteins(list(reported[protein_col]))
    #mapper = mapped_proteins[["primaryAccession", "uniProtKB_ID"]].set_index("primaryAccession").to_dict()["uniProtKB_ID"]
    #reported["protein"] = reported["Protein.Ids"].map(mapper)
    #return reported
    return mapped_proteins

def filter_and_map_triqler(triqler, fdr_threshold, fc_threshold, triqler_mapper):
    triqler_filtered = triqler[(triqler["q_value"] < fdr_threshold) & (abs(triqler["log2_fold_change"]) > fc_threshold )]
#    if triqler_mapper == False:
#        triqler_mapper = get_mapper(triqler_filtered, protein_col = "protein")
    triqler_filtered["uniProtKB_ID"] = triqler_filtered["protein"]
    triqler_filtered["geneName"] = triqler_filtered.uniProtKB_ID.map(triqler_mapper[["uniProtKB_ID", "geneName"]].set_index("uniProtKB_ID").to_dict()["geneName"])
    return triqler_filtered

def filter_and_map_reported_s3(reported_S3, fdr_threshold, fc_threshold, s3_mapper):
    reported_S3_filtered = reported_S3[(reported_S3["Adjusted_P_value"] < fdr_threshold) & (abs(reported_S3['log2（FC）']) > fc_threshold)]
#    if s3_mapper != False:
#        s3_mapper = get_mapper(reported_S3_filtered, protein_col = "Protein.Ids")
    reported_S3_filtered["uniProtKB_ID"] = reported_S3_filtered["Protein.Ids"].map(s3_mapper[["primaryAccession", "uniProtKB_ID"]].set_index("primaryAccession").to_dict()["uniProtKB_ID"])
    reported_S3_filtered["geneName"] = reported_S3_filtered.uniProtKB_ID.map(s3_mapper[["uniProtKB_ID", "geneName"]].set_index("uniProtKB_ID").to_dict()["geneName"])
    return reported_S3_filtered


def get_count_table(triqler, reported_S3):
    # Differential protein analysis
    overlapping_proteins = set(triqler["uniProtKB_ID"]) & set(reported_S3["uniProtKB_ID"]) 
    
    n_diff_triqler = len(triqler)
    n_diff_reported = len(reported_S3)
    n_diff_intersection = len(overlapping_proteins)
    n_diff_triqler_intersection = n_diff_triqler - n_diff_intersection
    n_diff_reported_intersection = n_diff_reported - n_diff_intersection
    n_diff_cols = ["triqler", "reported", "intersection", "triqler-intersection", "reported-intersection"]
    n_diff_protein_list = [n_diff_triqler, n_diff_reported, n_diff_intersection, n_diff_triqler_intersection, n_diff_reported_intersection]
    
    # remove na genes
    triqler = triqler[~triqler["geneName"].isna()]
    reported_S3 = reported_S3[~reported_S3["geneName"].isna()]
    
    overlapping_genes = set(triqler["geneName"]) & set(reported_S3["geneName"]) 
    n_diff_triqler = len(triqler)
    n_diff_reported = len(reported_S3)
    n_diff_intersection = len(overlapping_genes)
    n_diff_triqler_intersection = n_diff_triqler - n_diff_intersection
    n_diff_reported_intersection = n_diff_reported - n_diff_intersection
    
    n_diff_gene_list = [n_diff_triqler, n_diff_reported, n_diff_intersection, n_diff_triqler_intersection, n_diff_reported_intersection]
    
    count_table = pd.DataFrame([n_diff_protein_list, n_diff_gene_list], columns = n_diff_cols, index = ["diff_proteins", "diff_genes"])
    return count_table    

def enr_pathway_analysis(triqler, reported_S3, triqler_mapper, s3_mapper):
    triqler = triqler[~triqler["geneName"].isna()]
    reported_S3 = reported_S3[~reported_S3["geneName"].isna()]
    
    enr_reported = gp.enrichr(gene_list=reported_S3.geneName,
                     gene_sets=['KEGG_2019_Mouse'],
                     background=s3_mapper.geneName,
                     organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None, # don't write to disk
                    )
    enr_triqler = gp.enrichr(gene_list=triqler.geneName,
                     gene_sets=['KEGG_2019_Mouse'],
                     background=triqler_mapper.geneName,
                     organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None, # don't write to disk
                    )
    
    
    enr_reported.res2d["group"] = "reported"
    enr_triqler.res2d["group"] = "triqler"
    enr_clusters = pd.concat([enr_reported.res2d, enr_triqler.res2d])
    enr_clusters["pathway_DEG"] = enr_clusters.Overlap.map(lambda x:float(x.split("/")[0]))
    enr_clusters["pathway_genes"] = enr_clusters.Overlap.map(lambda x:float(x.split("/")[1]))
    enr_clusters["% Path is DEG"] = enr_clusters["pathway_DEG"] / enr_clusters["pathway_genes"]
    return enr_clusters


def add_pathways_to_count_table(count_table, enr_clusters):
    overlapping_pathways = intersection(enr_clusters[enr_clusters["group"] == "triqler"].Term.unique(), 
                 enr_clusters[enr_clusters["group"] == "reported"].Term.unique())
    n_pathways_triqler = len(enr_clusters[enr_clusters["group"] == "triqler"])
    n_pathways_reported = len(enr_clusters[enr_clusters["group"] == "reported"])
    n_pathways_overlapping = len(overlapping_pathways)
    n_pathways_triqler_intersection = n_pathways_triqler - n_pathways_overlapping
    n_pathways_reported_intersection = n_pathways_reported - n_pathways_overlapping
    n_pathways_list = [n_pathways_triqler, n_pathways_reported, n_pathways_overlapping, n_pathways_triqler_intersection, n_pathways_reported_intersection]

    count_table = count_table.T
    count_table["diff_pathways"] = n_pathways_list
    return count_table


def output_count_table(reported_S3, triqler, s3_mapper, triqler_mapper, fdr_threshold, fc_threshold,
                       pathway_fdr_threshold, percent_DEG_threshold):
    # S3_mapper is used as bgr for enrichr for S3 
    # triqler_mapper is used as bgr for enrichr for triqler
    
    #reported_S3 = pd.read_excel(reported_S3_file, header = 1)
    #triqler = parse_triqler(triqler_file)
    
    triqler = filter_and_map_triqler(triqler, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, triqler_mapper = triqler_mapper)
    reported_S3 = filter_and_map_reported_s3(reported_S3, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, s3_mapper = s3_mapper)

    count_table = get_count_table(triqler, reported_S3) # TABLE
    count_table_UP = get_count_table(triqler[triqler["log2_fold_change"]>0], reported_S3[reported_S3['log2（FC）']>0])
    count_table_DOWN = get_count_table(triqler[triqler["log2_fold_change"]<0], reported_S3[reported_S3['log2（FC）']<0])
    #count_table_down.rename(dict(zip(count_table_down.index,count_table_down.index.map(lambda x:x+"_Downregulated"))), axis = 1, inplace = True)

    enr_clusters = enr_pathway_analysis(triqler, reported_S3, triqler_mapper, s3_mapper)
    enr_clusters = enr_clusters[(enr_clusters["Adjusted P-value"] < pathway_fdr_threshold) &
                 (enr_clusters['% Path is DEG'] > percent_DEG_threshold)]


    count_table = add_pathways_to_count_table(count_table, enr_clusters)
    count_table.rename(dict(zip(count_table.columns,count_table.columns.map(lambda x:x+"_Total"))), axis = 1, inplace = True)    
 
    # upregulated
    enr_clusters_UP = enr_pathway_analysis(triqler[triqler["log2_fold_change"]>0], reported_S3[reported_S3['log2（FC）']>0], triqler_mapper, s3_mapper)
    enr_clusters_UP = enr_clusters_UP[(enr_clusters_UP["Adjusted P-value"] < pathway_fdr_threshold) &
                 (enr_clusters_UP['% Path is DEG'] > percent_DEG_threshold)]


    count_table_UP = add_pathways_to_count_table(count_table_UP, enr_clusters_UP)
    count_table_UP.rename(dict(zip(count_table_UP.columns,count_table_UP.columns.map(lambda x:x+"_Upregulated"))), axis = 1, inplace = True)

    # downregulated
    enr_clusters_DOWN = enr_pathway_analysis(triqler[triqler["log2_fold_change"]<0], reported_S3[reported_S3['log2（FC）']<0],triqler_mapper, s3_mapper)
    enr_clusters_DOWN = enr_clusters_DOWN[(enr_clusters_DOWN["Adjusted P-value"] < pathway_fdr_threshold) &
                 (enr_clusters_DOWN['% Path is DEG'] > percent_DEG_threshold)]

    count_table_DOWN = add_pathways_to_count_table(count_table_DOWN, enr_clusters_DOWN)
    count_table_DOWN.rename(dict(zip(count_table_DOWN.columns,count_table_DOWN.columns.map(lambda x:x+"_Downregulated"))), axis = 1, inplace = True)
    
    res = pd.concat([count_table, count_table_UP, count_table_DOWN], axis = 1)
    res["log2FC_fdr"] = fdr_threshold
    res["abs(log2FC)"] = fc_threshold
    res["pathway_fdr"] = pathway_fdr_threshold
    res["%DEG_in_pathway"] = percent_DEG_threshold
    return res
    #count_table.to_csv(output_name, sep = "\t")


def filter_enrichr(enr_clusters, pathway_fdr_threshold, percent_DEG_in_pathway_threshold):
    enr_clusters = enr_clusters[(enr_clusters["Adjusted P-value"] < pathway_fdr_threshold) & 
                                (enr_clusters['% Path is DEG'] > percent_DEG_in_pathway_threshold)]
    return enr_clusters
 
def check_inversed(triqler, reported_S3, fdr_threshold, fc_threshold, s3_mapper):

    reported_S3["uniProtKB_ID"] = reported_S3["Protein.Ids"].map(s3_mapper[["primaryAccession", "uniProtKB_ID"]].set_index("primaryAccession").to_dict()["uniProtKB_ID"])

    triqler_filtered = triqler[(triqler["q_value"] < fdr_threshold) & (abs(triqler["log2_fold_change"]) > fc_threshold )]
    reported_S3_filtered = reported_S3[(reported_S3["Adjusted_P_value"] < fdr_threshold) & (abs(reported_S3['log2（FC）']) > fc_threshold)]

    intersecting_proteins = intersection(triqler_filtered.protein,reported_S3_filtered["uniProtKB_ID"])

    triqler_filtered = triqler[triqler.protein.isin(intersecting_proteins)]
    reported_S3_filtered = reported_S3[reported_S3["uniProtKB_ID"].isin(intersecting_proteins)]
    
    triqler_filtered.set_index("protein", inplace = True)
    reported_S3_filtered.set_index("uniProtKB_ID", inplace = True)
    
    df_comparison = pd.DataFrame([triqler_filtered["log2_fold_change"], reported_S3_filtered['log2（FC）']]).T
    df_comparison.rename({"log2_fold_change":"triqler_FC", 'log2（FC）':"reported_FC"}, axis = 1, inplace = True)
    df_comparison["inversed"] = df_comparison["triqler_FC"] * df_comparison["reported_FC"] < 0
    return df_comparison


os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-09-22_ctrl_vs_ST_vs_LT_study")


def read_CtrlST(reported_file = "S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx"):
    return pd.read_excel(reported_file, header = 1).rename({"ST-Ctrl_diff":'log2（FC）',"ST-Ctrl_p adj":"Adjusted_P_value"}, axis = 1)
    

def read_CtrlLT(reported_file = "S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx"):
    return pd.read_excel(reported_file, header = 1).rename({"LT-Ctrl_diff":'log2（FC）',"LT-Ctrl_p adj":"Adjusted_P_value"}, axis = 1)


def read_LTST(reported_file = "S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx"):
    return pd.read_excel(reported_file, header = 1).rename({"LT-ST_diff":'log2（FC）',"LT-ST_p adj":"Adjusted_P_value"}, axis = 1)


##############################
# Get protein specific table #
##############################

fdr_threshold = 1
fc_threshold = -5
pathway_fdr_threshold = 0.05
percent_DEG_in_pathway_threshold = 0.0

# LTST 
reported = read_LTST() #########
triqler = parse_triqler("LT_ST/proteins_fc_0.1") #########

triqler[triqler.protein == "CASP3_MOUSE"].log2_fold_change


s3_mapper = get_mapper(reported, protein_col = "Protein.Ids")
triqler_mapper = get_mapper(triqler, protein_col = "protein")


#triqler["log2_fold_change"] = -triqler["log2_fold_change"]
check_inversed(triqler, reported, fdr_threshold, fc_threshold, s3_mapper) # just a check

triqler = filter_and_map_triqler(triqler, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, triqler_mapper = triqler_mapper)
reported = filter_and_map_reported_s3(reported, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, s3_mapper = s3_mapper)
reported = reported.dropna(subset=['uniProtKB_ID'])

triqler_FC = triqler.set_index("uniProtKB_ID")[["log2_fold_change", "q_value"]].rename({"log2_fold_change":"triqler_log2FC", "q_value":"triqler_FDR"}, axis =1)
reported_FC = reported.set_index("uniProtKB_ID")[['log2（FC）', 'Adjusted_P_value']].rename({'log2（FC）':"reported_log2FC", "Adjusted_P_value": "reported_FDR"}, axis = 1)


res = pd.concat([triqler_FC, reported_FC], axis = 1)
res = res.rename(dict(zip(res.columns,[i+"_LT_ST" for i in res.columns])), axis = 1)

LTST = res 
res.to_csv("LT_ST.csv", sep = "\t") #########

# CTRLST

reported = read_CtrlST() #########
triqler = parse_triqler("ctrl_ST/proteins_fc_0.1") #########

triqler["log2_fold_change"] = -triqler["log2_fold_change"]
check_inversed(triqler, reported, fdr_threshold, fc_threshold, s3_mapper) # just a check

triqler = filter_and_map_triqler(triqler, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, triqler_mapper = triqler_mapper)
reported = filter_and_map_reported_s3(reported, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, s3_mapper = s3_mapper)
reported = reported.dropna(subset=['uniProtKB_ID'])

triqler_FC = triqler.set_index("uniProtKB_ID")[["log2_fold_change", "q_value"]].rename({"log2_fold_change":"triqler_log2FC", "q_value":"triqler_FDR"}, axis =1)
reported_FC = reported.set_index("uniProtKB_ID")[['log2（FC）', 'Adjusted_P_value']].rename({'log2（FC）':"reported_log2FC", "Adjusted_P_value": "reported_FDR"}, axis = 1)


res = pd.concat([triqler_FC, reported_FC], axis = 1)
res = res.rename(dict(zip(res.columns,[i+"_CTRL_ST" for i in res.columns])), axis = 1)
CTRLST = res
res.to_csv("Ctrl_ST.csv", sep = "\t") #########


# CTRLLT

reported = read_CtrlLT() #########
triqler = parse_triqler("ctrl_LT/proteins_fc_0.1") #########

triqler["log2_fold_change"] = -triqler["log2_fold_change"]
check_inversed(triqler, reported, fdr_threshold, fc_threshold, s3_mapper) # just a check

triqler = filter_and_map_triqler(triqler, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, triqler_mapper = triqler_mapper)
reported = filter_and_map_reported_s3(reported, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, s3_mapper = s3_mapper)
reported = reported.dropna(subset=['uniProtKB_ID'])

triqler_FC = triqler.set_index("uniProtKB_ID")[["log2_fold_change", "q_value"]].rename({"log2_fold_change":"triqler_log2FC", "q_value":"triqler_FDR"}, axis =1)
reported_FC = reported.set_index("uniProtKB_ID")[['log2（FC）', 'Adjusted_P_value']].rename({'log2（FC）':"reported_log2FC", "Adjusted_P_value": "reported_FDR"}, axis = 1)


res = pd.concat([triqler_FC, reported_FC], axis = 1)
res = res.rename(dict(zip(res.columns,[i+"_CTRL_LT" for i in res.columns])), axis = 1)
CTRLLT = res
res.to_csv("Ctrl_LT.csv", sep = "\t") #########



final = pd.concat([CTRLST, CTRLLT, LTST], axis = 1)
final = final[~final.index.str.contains("DECOY")]
final.to_csv("protein_table.tsv", sep = "\t")



uniProt_to_KEGG_mapper = get_uniProtKB_ID_to_KEGG_mapper(final.index)
final["KEGG_id"] = final.index.map(uniProt_to_KEGG_mapper)

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST


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

# pathway related proteins...
    
pathways = oxaliplatin_pathways()

pathways
apoptosis = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "Apoptosis")
apoptosis.to_csv("apoptosis.tsv", sep = "\t")
nucleotide_excision_repair = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "Nucleotide excision repair")
nucleotide_excision_repair.to_csv("nucleotide_excision_repair.tsv", sep = "\t")
mismatch_repair = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "Mismatch repair")
mismatch_repair.to_csv("mismatch_repair.tsv", sep = "\t")
erbB_signaling_pathway = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "ErbB signaling pathway")
erbB_signaling_pathway.to_csv("erbB_signaling_pathway.tsv", sep = "\t")
cell_cycle = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "Cell cycle")
cell_cycle.to_csv("cell_cycle.tsv", sep = "\t")
p53 = get_KEGG_pathway_proteins(final_df = final, pathway_df = pathways, pathway_term = "p53 signaling pathway")
p53.to_csv("p53_signaling_pathway.tsv", sep = "\t")
 

final.to_csv("protein_table_output.tsv", sep = "\t")













