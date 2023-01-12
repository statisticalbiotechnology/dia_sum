#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 18:50:31 2022

@author: ptruong
"""
import os
import pandas as pd
from time import time
import argparse



parser = argparse.ArgumentParser(
    description='Concatenates osw output to one single file.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--file_directory', type=str,
                    help='OSW results file directory.')

parser.add_argument('--output', type=str, default="concatenated_osw_results.csv", 
                    help='output name.')

# parse arguments from command line
args = parser.parse_args()
file_directory = args.file_directory
output = args.output

def concatenate_osw_results(file_directory, output):
    start = time()
    dfs = []
    for i in os.listdir(file_directory):
        if i[-10:] == "dscore.csv":
            df = pd.read_csv(file_directory + "/" +  i, sep = "\t")
            dfs.append(df)
            print(time()-start)
    end = time()
    print(end-start)
    
    df = pd.concat(dfs)
    df.Intensity = df.Intensity.map(lambda x: x if (type(x) == float) else float(x.replace(",", ".")))
    df.ProteinName = df.ProteinName.map(lambda x:x.replace("YEAS8", "YEAST"))
    df["ProteinNameCount"] = df.ProteinName.map(lambda x:len(x.split("|"))) # removed PSMS that have more than 1 protein
    df = df[df.ProteinNameCount <= 3]
    df.to_csv(output, sep = "\t")

if __name__ == "__main__":
    concatenate_osw_results(file_directory, output)
    
    
