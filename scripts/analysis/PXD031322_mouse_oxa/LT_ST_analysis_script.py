#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 19:19:18 2022

@author: ptruong
"""

import pandas as pd 
import numpy as np
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *
import gseapy as gp


paper_term = ["Ribosome", "Spliceosome", "Endocytosis", "Steroid biosynthesis",
              "Dopaminergic synapse", "RNA transport", "Glutathione metabolism",
              "Proteasome", "Tight junction", "Complement and coagulation cascades", 
              "Metabolism of xenobiotics by cytochrome P450", "Fructose and mannose metabolism",
              "mRNA surveillance pathway", "Arginine biosynthesis", "Protein export",
              "Protein processing in endoplasmic reticulum", "N-Glycan biosynthesis",
              "Tyrosine metabolism", "Metabolic pathways", "Adrenaergic signaling in cardiomyocytes"]

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


os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-09-22_ST_vs_LT_study")

# divide into S3 and S4 analysis and into 0.1 to 0.6 FC, including fc 0.415 
# Count protein [Done]
# Protein overlap [Done]
# Perform pathway analysis []

fc_threshold = 0.2
fdr_threshold = 0.05

reported_S3 = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
reported_S3.rename({"LT-ST_p adj":"Adjusted_P_value", "LT-ST_diff":"log2（FC）"}, axis = 1, inplace = True)
#reported_S4 = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
#reported_S5 = pd.read_excel("S5_KEGG_enrichment_analysis_of_proteins_in_six_subclusters.xlsx", header = 1) # thier Kegg pathway, as reference
triqler = parse_triqler("proteins_fc_0.1")
s3_mapper = get_mapper(reported_S3, protein_col = "Protein.Ids")
triqler_mapper = get_mapper(triqler, protein_col = "protein")
#S3_statistics(reported_S3)
#S4_statistics(reported_S4)
#os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-09-20_ctrl_vs_LT_study")

fc_threshold = 0.1
fdr_threshold = 0.05

triqler = parse_triqler("proteins_fc_0.1")
output_name = f"total_count_table_fc_{fc_threshold}_fdr_{fdr_threshold}.tsv"
reported_S3 = pd.read_excel("S3_list_of_proteins_differentially_regulated_between_LT_and_Ctrl.xlsx", header = 1)
#reported_S4 = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
#reported_S5 = pd.read_excel("S5_KEGG_enrichment_analysis_of_proteins_in_six_subclusters.xlsx", header = 1) # thier Kegg pathway, as reference

fdr_threshold = 0.05
fc_threshold = 0.1
pathway_fdr_threshold = 0.05
percent_DEG_in_pathway_threshold = 0.0
reported_S3 = pd.read_excel("S4_list_of_differential_regulation_protein_in_six_subclusters.xlsx", header = 1)
reported_S3.rename({"LT-ST_p adj":"Adjusted_P_value", "LT-ST_diff":"log2（FC）"}, axis = 1, inplace = True)
triqler = parse_triqler("proteins_fc_0.1")
triqler["log2_fold_change"] = -triqler["log2_fold_change"]
#check = check_inversed(triqler, reported_S3, fdr_threshold, fc_threshold, s3_mapper)
res = output_count_table(reported_S3, triqler, s3_mapper = s3_mapper, triqler_mapper = triqler_mapper, 
           fdr_threshold = fdr_threshold, fc_threshold = fc_threshold,
           pathway_fdr_threshold = pathway_fdr_threshold, 
           percent_DEG_threshold = percent_DEG_in_pathway_threshold)
res.to_csv(f"count_table_fc_{fc_threshold}_fdr_{fdr_threshold}_pathwayFDR_{pathway_fdr_threshold}_percentDEG_{percent_DEG_in_pathway_threshold}.tsv", sep = "\t")


#alysis - biological

fdr_threshold = 0.05
fc_threshold = 0.1
triqler = filter_and_map_triqler(triqler, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, triqler_mapper = triqler_mapper)
reported_S3 = filter_and_map_reported_s3(reported_S3, fdr_threshold = fdr_threshold, fc_threshold = fc_threshold, s3_mapper = s3_mapper)
enr_clusters = enr_pathway_analysis(triqler, reported_S3, triqler_mapper, s3_mapper)

enr_triqler = enr_clusters[enr_clusters.group == "triqler"]
enr_triqler["triqler_pathway_FDR"] = enr_triqler["Adjusted P-value"]
enr_reported = enr_clusters[enr_clusters.group == "reported"]
enr_reported["reported_pathway_FDR"] = enr_reported["Adjusted P-value"]
enr_triqler.set_index("Term", inplace = True)
enr_reported.set_index("Term", inplace = True)
comparison = pd.concat([enr_triqler["triqler_pathway_FDR"], enr_reported["reported_pathway_FDR"]], axis = 1)

#comparison = comparison.where(comparison < 0.05) #filter away all values with higher than 0.05 fdr
comparison = comparison.dropna(subset=["triqler_pathway_FDR", "reported_pathway_FDR"], how  = "all")
comparison = comparison.replace(np.nan, 1)

comparison["more_significant"] = (comparison["triqler_pathway_FDR"] - comparison["reported_pathway_FDR"]) > 0

def method_is_significant(x):
    if x == False:
        return "triqler"
    else:
        return "reported"

comparison["more_significant"] = comparison["more_significant"].map(lambda x:method_is_significant(x))
comparison = comparison.replace(1,np.nan)
# FIX THIS LIST!!
# Remove nan columns
# and fix the signal column
comparison[comparison.more_significant == "triqler"]
comparison[comparison.more_significant == "reported"]

comparison[~comparison.triqler_pathway_FDR.isna()]
comparison[~comparison.reported_pathway_FDR.isna()]

comparison.to_csv(f"pathway_comparison_fcFDR_{fdr_threshold}_fcThreshold_{fc_threshold}.tsv", sep = "\t")
comparison.to_csv(f"pathway_comparison_fcFDR_{fdr_threshold}_fcThreshold_{fc_threshold}_noFilter.tsv", sep = "\t")

