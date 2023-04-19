#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:09:35 2022

@author: ptruong

https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html

"""

# Standard library packages
import io
import os
import pandas as pd

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
#from Bio.KEGG.KGML import KGML_parser
#from Bio.Graphics.KGML_vis import KGMLCanvas




# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# Some code to return a Pandas dataframe, given tabular text
def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

# Find a specific entry with a precise search term
#result = REST.kegg_get("path:mmu04210").read()
#print(result)

# HERE WE HAVE INFORMATION!!!!
#mmu04210

#https://rest.kegg.jp/link/mmu/mmu04210


# here we have genes related to apoptosis

def get_mouse_pathway(pathway_id = "path:mmu04210"):
    result = REST.kegg_link(target_db = "mmu", source_db = pathway_id).read()
    res = to_df(result)
    res.rename({0:"pathway", 1:"gene"}, axis = 1)
    return res
