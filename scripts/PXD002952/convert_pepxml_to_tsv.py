#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 00:27:03 2021

@author: ptruong
"""

import os 

from pyteomics import pepxml

os.chdir("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger") 

pxml = pepxml.DataFrame(file)

for file in os.listdir():
    if file[:8] == 'interact':
        if file[-4:] == ".xml":
            print(file)
            pxml = pepxml.DataFrame(file)
            pxml.to_csv(file+".tsv", sep = "\t", index = False)


















