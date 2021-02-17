#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 22:28:28 2021

@author: ptruong
"""

import pandas as pd
import numpy as np

def gen_swath_windows(start = 350, end = 1400, window_size = 50):
    windows = []
    for i in range(start,end+1, window_size):
        windows.append(i)
    
    
    windows = np.array(windows)
    
    
    swath_window_df = pd.DataFrame(np.array([windows[:-1], windows[1:]]).T, columns = ["start", "end"])
    return swath_window_df


swath_window_df = gen_swath_windows(start = 350, end = 1400, window_size = 49)

swath_window_df.to_csv("swath_window_49.csv", sep = "\t", index = False)
