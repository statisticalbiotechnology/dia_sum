#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:41:49 2022

@author: ptruong
"""

import pandas as pd 
import argparse


def copy_top3_to_scatter_input(input_file, output):
    df = pd.read_csv(input_file, sep = "\t")
    df["log2(A,B)"] = -df["log2(A,B)"]
    df.to_csv(output, sep="\t", index = False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='File copies top3_results file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--input_file', type=str,
                        help='Scatterplot converter input file.')

    parser.add_argument('--output', type=str,
                        help='Output name.')


    # parse arguments from command line
    args = parser.parse_args()
    input_file = args.input_file
    output = args.output
    
    copy_top3_to_scatter_input(input_file, output)
    
    