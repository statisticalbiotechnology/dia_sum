#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 02:23:48 2020

@author: ptruong
"""


def count_melted_for_all_samples(df, specie = None):
    """
    df = melted triq or spec.
    """
    samples = ["S0"+str(i) for i in range(1,10)] + ["S10"]
    counts = []
    for i in samples:
        count_df = df[df["sample"] == i]
        if specie == "ARATH":
            count_df = count_df[count_df["specie"] == "ARATH"].dropna()
        elif specie == "HUMAN":
            count_df = count_df[count_df["specie"] == "HUMAN"].dropna()
        elif specie == "CAEEL":
            count_df = count_df[count_df["specie"] == "CAEEL"].dropna()
        count = len(count_df.dropna())
        counts.append(count)
    return counts

