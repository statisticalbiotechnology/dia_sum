#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 12:16:23 2022

@author: ptruong
"""

import pandas as pd
import numpy as np


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

def get_protein_mapper(fasta = "swissprot_mouse_20210115.fasta"):
    filename = fasta
    protein_name_1_list = []
    protein_name_2_list = []
    with open(filename, 'r', encoding='UTF-8') as file:
        while (line := file.readline().rstrip()):
            if line[0] == ">":
                protein_name_1 = line.split("|")[1]
                protein_name_2 = line.split("|")[2].split(" ")[0]
                protein_name_1_list.append(protein_name_1)
                protein_name_2_list.append(protein_name_2)
    protein_map = dict(zip(protein_name_1_list, protein_name_2_list))
    protein_map_reversed = dict(zip(protein_name_2_list, protein_name_1_list))
    return protein_map, protein_map_reversed

# fix protein groups
def get_triqler_protein_groups(triqler):
    protein_groups_triqler = triqler[triqler.protein.str.contains(";")]
    pg_list = []
    for pg in protein_groups_triqler.protein:
        for protein in pg.split(";"):
            if protein not in pg_list:
                pg_list.append(protein)
    
    proteins_that_are_uniquely_identified_in_pg = triqler[triqler.protein.isin(pg_list)].protein
    
    pg_substr_list= [] #If any protein och protein group constains anything from protein accession without isoform number, we remove it from the complete list.
    for protein in proteins_that_are_uniquely_identified_in_pg:
        protein_substr = protein.split("_")[0][:-1]
        if protein_substr not in pg_substr_list:
            pg_substr_list.append(protein_substr)
    
    pg_keep_list = []
    for pg in protein_groups_triqler.protein:
        for protein in pg.split(";"):
            protein_substr = protein.split("_")[0][:-1]
            if protein_substr not in pg_substr_list:
                if pg not in pg_keep_list:
                    pg_keep_list.append(pg)
    return protein_groups_triqler, pg_list, proteins_that_are_uniquely_identified_in_pg, pg_substr_list, pg_keep_list


# Do Venn-diagram between proteins in triqler and reported in article 
# Triqler has 279 differentially abundant with fc = 0.4
# Reported har 280 differentially abundant at fc = 0.58

#triqler_results_dir = "2022-08-11_run/2022-08-22_triqler_inputs/proteins_fc_0.415_triqler_input_report_ST_CTRL.tsv"
#triqler_results_dir = "2022-08-11_run/2022-08-22_triqler_inputs/proteins_fc_0.415_triqler_input_report_LT_CTRL.tsv"

triqler_results_dir = "2022-08-11_run/2022-08-22_triqler_inputs/proteins_fc_0.415_triqler_input_report_ST_LT_CTRL.1vs3.CTRL_vs_ST.tsv"
#triqler_results_dir = "2022-08-11_run/2022-08-22_triqler_inputs/proteins_fc_0.415_triqler_input_report_ST_LT_CTRL.1vs2.CTRL_vs_LT.tsv"

report_results = "reported_ST_vs_CTRL.xlsx"
#report_results = "reported_LT_vs_CTRL.xlsx"
fasta_file = "swissprot_mouse_20210115.fasta"
fdr_threshold = 0.05

triqler = parse_triqler(triqler_results_dir)
triqler = triqler[triqler.q_value < fdr_threshold]
reported = pd.read_excel(report_results, skiprows=[0])
reported = reported[~reported.Condition.isin(["Not_sig"])]
protein_map, protein_map_reversed = get_protein_mapper(fasta = fasta_file)
reported["protein"] = reported["Protein.Ids"].map(protein_map)
reported = reported[reported["Adjusted_P_value"] < fdr_threshold]
protein_groups_triqler, pg_list, proteins_that_are_uniquely_identified_in_pg, pg_substr_list, pg_keep_list = get_triqler_protein_groups(triqler)



# Remove PG: who has uniquely identified proteoforms (e.g. without proteoform number information)
# E.g. if we have MYH2, we remove MYH1;MYH2; MYH8 etc.
triqler_keep_pg = triqler[(triqler.protein.str.contains(";") & triqler.protein.isin(pg_keep_list))]
triqler_non_pg = triqler[~triqler.protein.str.contains(";")]
triqler = pd.concat([triqler_non_pg, triqler_keep_pg])

pg_proteins = []
for pg in triqler_keep_pg.protein:
    for protein in pg.split(";"):
        pg_proteins.append(protein)
number_of_overlap_in_pg_proteins = len(reported[reported.protein.isin(pg_proteins)])
print(number_of_overlap_in_pg_proteins) #0 we can ignore this part of the analysis.



overlapping_differential_abundant_proteins = len((set(triqler.protein) & set(reported.protein)))
n_reported_protein_non_overlap = len(reported.protein) - overlapping_differential_abundant_proteins
n_triqler_protein_non_overlap = len(triqler.protein) - overlapping_differential_abundant_proteins
n_reported = len(reported.protein)
n_triqler = len(triqler.protein)

# library
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

overlapping_differential_abundant_proteins
n_reported_protein_non_overlap
n_triqler_protein_non_overlap
n_reported
n_triqler

venn2(subsets = (n_triqler_protein_non_overlap,
                 n_reported_protein_non_overlap,
      overlapping_differential_abundant_proteins),
      set_labels = (f'Triqler ({n_triqler} proteins)', f'Reported ({n_reported} proteins)'))
plt.title(f"(LT vs CTRL) Venn diagram of differentially abundant proteins, FDR < {fdr_threshold}, fc = 0.415")
plt.show()


#### varying differential expression...

fdr_threshold = 0.05
triqler = parse_triqler(triqler_results_dir)
triqler = triqler[triqler.protein_id_posterior_error_prob < fdr_threshold]
reported = pd.read_excel(report_results, skiprows=[0])
reported = reported[~reported.Condition.isin(["Not_sig"])]
protein_map, protein_map_reversed = get_protein_mapper(fasta = fasta_file)
reported["protein"] = reported["Protein.Ids"].map(protein_map)
reported = reported[reported["Adjusted_P_value"] < fdr_threshold]
protein_groups_triqler, pg_list, proteins_that_are_uniquely_identified_in_pg, pg_substr_list, pg_keep_list = get_triqler_protein_groups(triqler)

reported.columns
reported["abs(log2FC)"] = abs(reported['log2（FC）'])
reported_fc = reported[(reported['abs(log2FC)'] > 0.4)]
reported_fc[reported_fc["Adjusted_P_value"] < 0.05]
reported_fc[reported_fc["Condition"] != "Not_sig"]

# These parameters are just not right?
reported[reported['log2（FC）']>1.5]
reported[(reported['log2（FC）']<0.67) & (reported["Adjusted_P_value"]<0.05) ] # THIS IS NOT CORRECT? 

#################






















