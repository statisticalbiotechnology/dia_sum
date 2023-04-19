#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:48:56 2022

@author: ptruong
"""

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


df["condition"] = df.Run.map(lambda x:x.split("-"[-1])[-1][:-1])
df.columns



# Cluster on triqlers quantities


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

def quant_norm(df):
    ranks = (df.rank(method="first")
              .stack())
    rank_mean = (df.stack()
                   .groupby(ranks)
                   .mean())
    # Add interpolated values in between ranks
    finer_ranks = ((rank_mean.index+0.5).to_list() +
                    rank_mean.index.to_list())
    rank_mean = rank_mean.reindex(finer_ranks).sort_index().interpolate()
    return (df.rank(method='average')
              .stack()
              .map(rank_mean)
              .unstack())
#quant_norm(df)

df = parse_triqler("proteins.1vs2.tsv")
cols = ['peptides','q_value','posterior_error_prob','num_peptides',
        'protein_id_posterior_error_prob','log2_fold_change', 'diff_exp_prob_0.6']
vals = df.iloc[:,~df.columns.isin(cols)].set_index("protein")



import seaborn as sns

sns.heatmap(vals)
df.drop("protein")
len(clusters.labels_)

vals.groupby("cluster").count()

vals
vals_clr.max()
vals_clr.min()

vals.replace(1,np.nan)

### Clustering on relative abundances

import numpy as np
from skbio.stats.composition import clr

vals

x = np.array([.1, .3, .4, .2])
vals_clr = clr(vals)
vals_clr = pd.DataFrame(vals_clr, index = vals.index, columns = vals.columns)

vals.replace(1,np.nan)

from sklearn.cluster import AgglomerativeClustering
clusters = AgglomerativeClustering(n_clusters=6, affinity='euclidean', linkage='complete')
X = vals_clr.values
clusters.fit_predict(X)
vals_clr["cluster"] = clusters.labels_
vals_clr.groupby("cluster").count()
vals

####### THIS DID NOT WORK...







