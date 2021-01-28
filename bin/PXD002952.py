#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 20:19:50 2021

@author: ptruong
"""
import os 

import pandas as pd 


files = os.listdir()


df = pd.read_csv("Spectronaut_onlyHumanPeptidesDetectedByDIAUmpire_Report_ProtHUMAN.tsv", sep = "\t")


df.columns

run = df[""]
df_triqler_formatted = pd.DataFrame()












