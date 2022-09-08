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
import gseapy as gp
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *
from gseapy_plot import *

os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322/2022-08-11_run")


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

def qprofiler_run(query_genes, background_genes, user_threshold = 0.05):
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
        'background':background_genes
        },
        headers={
        'User-Agent':'FullPythonRequest'
        }
    )
    
    res = r.json()['result']
    return res

def plot_dotplot(df, size = 10, title = "KEGG"):
    temp = df["Overlap"].str.split("/", expand=True).astype(int)
    df = df.assign(Hits_ratio=temp.iloc[:, 0] / temp.iloc[:, 1])
    # make area bigger to better visualization
    # area = df["Hits_ratio"] * plt.rcParams["lines.linewidth"] * 100
    area = np.pi * (df["Hits_ratio"] * size * plt.rcParams["lines.linewidth"]).pow(2)
    
    xlabel = "Cluster"
    
    x = df["group"]
    ylabels = df["Term"].values
    cbar_title = r"Adjusted P-value"
    figsize: Tuple[float] = (6, 5.5)
    fig, ax = plt.subplots(figsize=figsize)
    colname = "Adjusted P-value"
    #colmap = df[colname].round().astype("int")
    colmap = df[colname]
    vmin = np.percentile(colmap.min(), 2)
    vmax = np.percentile(colmap.max(), 98)
    #cmap: str = "viridis_r"
    cmap: str = "viridis"
    sc = ax.scatter(
        x=x,
        y=ylabels,
        s=area,
        edgecolors="face",
        c=colmap,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    
    
    ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
    ax.xaxis.set_tick_params(labelsize=14)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.set_axisbelow(True)  # set grid blew other element
    ax.grid(axis="y")  # zorder=-1.0
    ax.margins(x=0.25)
    
    # We change the fontsize of minor ticks label
    # ax.tick_params(axis='y', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=14)
    
    # scatter size legend
    # we use the *func* argument to supply the inverse of the function
    # used to calculate the sizes from above. The *fmt* ensures to string you want
    handles, labels = sc.legend_elements(
        prop="sizes",
        num=3,  # fmt="$ {x:.2f}",
        color="gray",
        func=lambda s: np.sqrt(s / np.pi) / plt.rcParams["lines.linewidth"] / size,
    )
    ax.legend(
        handles,
        labels,
        title="% Path\nis DEG",
        bbox_to_anchor=(1.02, 0.9),
        loc="upper left",
        frameon=False,
    )
    # colorbar
    # cax = fig.add_axes([1.0, 0.20, 0.03, 0.22])
    cbar = fig.colorbar(
        sc,
        shrink=0.2,
        aspect=10,
        anchor=(0.0, 0.2),  # (0.0, 0.2),
        location="right"
        # cax=cax,
    )
    # cbar.ax.tick_params(right=True)
    cbar.ax.set_title(cbar_title, loc="left", fontweight="bold")
    for key, spine in cbar.ax.spines.items():
        spine.set_visible(False)
    
    ax.set_title(title, fontsize=20, fontweight="bold")

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

bgr_mapped = get_mapped_proteins(ids = list(bgr_proteins))
# we might need to chunk this for api request quicker



"""
user_threshold = 0.05

c1_pathways = qprofiler_run(query_genes = list(c1_mapped.geneName),
              background_genes = list(c1_bgr_mapped.geneName),
              user_threshold = 0.05)
c2_pathways = qprofiler_run(query_genes = list(c2_mapped.geneName),
              background_genes = list(c2_bgr_mapped.geneName),
              user_threshold = 0.05)
c3_pathways = qprofiler_run(query_genes = list(c3_mapped.geneName),
              background_genes = list(c3_bgr_mapped.geneName),
              user_threshold = 0.05)
c4_pathways = qprofiler_run(query_genes = list(c4_mapped.geneName),
              background_genes = list(c4_bgr_mapped.geneName),
              user_threshold = 0.05)
c5_pathways = qprofiler_run(query_genes = list(c5_mapped.geneName),
              background_genes = list(c5_bgr_mapped.geneName),
              user_threshold = 0.05)
c6_pathways = qprofiler_run(query_genes = list(c6_mapped.geneName),
              background_genes = list(c6_bgr_mapped.geneName),
              user_threshold = 0.05)

"""


#gene_list = c1_mapped.geneName
#background_list = c1_bgr_mapped.geneName
mouse_library = gp.get_library_name(organism="Mouse")

enr_c1 = gp.enrichr(gene_list=c1_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
enr_c2 = gp.enrichr(gene_list=c2_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
enr_c3 = gp.enrichr(gene_list=c3_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
enr_c4 = gp.enrichr(gene_list=c4_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
enr_c5 = gp.enrichr(gene_list=c5_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )
enr_c6 = gp.enrichr(gene_list=c6_mapped.geneName,
                 gene_sets=['KEGG_2019_Mouse'],
                 background=bgr_mapped.geneName,
                 organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

enr_c1.res2d["group"] = "C1"
enr_c2.res2d["group"] = "C2"
enr_c3.res2d["group"] = "C3"
enr_c4.res2d["group"] = "C4"
enr_c5.res2d["group"] = "C5"
enr_c6.res2d["group"] = "C6"

enr_clusters = pd.concat([enr_c1.res2d, enr_c2.res2d, enr_c3.res2d,enr_c4.res2d, enr_c5.res2d, enr_c6.res2d])
enr_clusters["pathway_DEG"] = enr_clusters.Overlap.map(lambda x:float(x.split("/")[0]))
enr_clusters["pathway_genes"] = enr_clusters.Overlap.map(lambda x:float(x.split("/")[1]))
enr_clusters["% Path is DEG"] = enr_clusters["pathway_DEG"] / enr_clusters["pathway_genes"]

paper_term = ["Ribosome", "Spliceosome", "Endocytosis", "Steroid biosynthesis",
              "Dopaminergic synapse", "RNA transport", "Glutathione metabolism",
              "Proteasome", "Tight junction", "Complement and coagulation cascades", 
              "Metabolism of xenobiotics by cytochrome P450", "Fructose and mannose metabolism",
              "mRNA surveillance pathway", "Arginine biosynthesis", "Protein export",
              "Protein processing in endoplasmic reticulum", "N-Glycan biosynthesis",
              "Tyrosine metabolism", "Metabolic pathways", "Adrenaergic signaling in cardiomyocytes"]

plot_input = enr_clusters[enr_clusters["Adjusted P-value"] < 0.05]
plot_input_paper_term = plot_input[plot_input.Term.isin(paper_term)]

df = plot_input
size = 10 # dot size baseline
title="KEGG_2019_mouse" #title

plot_dotplot(df = plot_input, size = 15, title = "KEGG_2019_mouse, fc_eval = 0.415, fdr = 0.05, all pathways")
plot_dotplot(df = plot_input_paper_term, size = 15, title = "KEGG_2019_mouse, fc_eval = 0.415, fdr = 0.05, reported pathways")


df.to_csv("enrichr_triqler_fc_0.415_whole_bgr.tsv", sep = "\t", index = False)







