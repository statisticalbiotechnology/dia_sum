#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:00:38 2022

@author: ptruong
"""

import os


import pandas as pd 
import numpy as np
import gseapy as gp
import os
os.chdir("/home/ptruong/git/dia_sum/scripts/clean/PXD031322_mouse_oxa")
from uniprot_idmapper import *
from gseapy_plot import *

os.chdir("/hdd_14T/data/PXD031322_oxaliplatin_dia_study/ftp.pride.ebi.ac.uk/pride/data/archive/2022/07/PXD031322")


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

paper_term = ["Ribosome", "Spliceosome", "Endocytosis", "Steroid biosynthesis",
              "Dopaminergic synapse", "RNA transport", "Glutathione metabolism",
              "Proteasome", "Tight junction", "Complement and coagulation cascades", 
              "Metabolism of xenobiotics by cytochrome P450", "Fructose and mannose metabolism",
              "mRNA surveillance pathway", "Arginine biosynthesis", "Protein export",
              "Protein processing in endoplasmic reticulum", "N-Glycan biosynthesis",
              "Tyrosine metabolism", "Metabolic pathways", "Adrenaergic signaling in cardiomyocytes"]



df = pd.read_excel("reported_6_subclusters.xlsx", header = 1)
bgr_genes = df["Gene.Symbol"]

cluster_dict = {}
for cluster in df[df.Condition != "Others"].Condition.unique():
    
    df_cluster = df[df.Condition == cluster]
    
    print(cluster)
    print("ST-CTRL p-value max" + str(df_cluster["ST-Ctrl_p adj"].max()))
    print("LT-CTRL p-value max" + str(df_cluster["LT-Ctrl_p adj"].max()))
    print("LT-ST p-value max" + str(df_cluster["LT-ST_p adj"].max()))
    print("")
    
    cluster_dict[cluster] = df[df.Condition == cluster]



enr_cluster_list = []
for cluster in df[df.Condition != "Others"].Condition.unique():    
    enr_cluster = gp.enrichr(gene_list=cluster_dict[cluster]["Gene.Symbol"][~cluster_dict[cluster]["Gene.Symbol"].isna()],
                     gene_sets=['KEGG_2019_Mouse'],
                     background=bgr_genes,
                     organism='mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None, # don't write to disk
                    )
    enr_cluster.res2d["group"] = cluster
    enr_cluster_list.append(enr_cluster.res2d)




enr_clusters = pd.concat(enr_cluster_list)
fdr_threshold = 0.05
plot_input = enr_clusters[enr_clusters["Adjusted P-value"] < 0.05]
plot_input_paper_term = plot_input[plot_input.Term.isin(paper_term)]

df = plot_input
size = 10 # dot size baseline
title="KEGG_2019_mouse" #title

plot_dotplot(df = plot_input, size = 15, title = "KEGG_2019_mouse, fdr = 0.05, all pathways")
plot_dotplot(df = plot_input_paper_term, size = 15, title = "KEGG_2019_mouse, fdr = 0.05. reported pathways")

# Missing relative the paper
# ["Metabolic pathways", "Adrenaergic signaling in cardiomyocytes"]
df.to_csv("enrichr_reported.tsv", sep = "\t", index = False)


