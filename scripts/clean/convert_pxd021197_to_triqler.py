#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:07:05 2022

@author: ptruong
"""

import pandas as pd
import argparse

def control_condition_rename(x):
    if x in ["Healthy control 1", "Healthy control 2", "Healthy control 3", "Healthy control 4", "Bacterial pneumonia control"]:
        return "control"
    else:
        return x
    
def convert_pxd021197_to_triqler(filename = "plasma.aligned.tsv", mapper_file = "PXD021197_DIA_DDA.xlsx", output = "output.tsv"):
    mapper = pd.read_excel(mapper_file)
    df = pd.read_csv(filename, sep = "\t")
    df["mapper"] = df.filename.map(lambda x:x.split(".")[0])
    mapper.SampleName
    len(mapper.Patient.unique())
    mapper.columns
    patient_map=mapper[["FileName", "Patient"]].set_index("FileName").to_dict()["Patient"]
    run_map=mapper[["FileName", "SampleName"]].set_index("FileName").to_dict()["SampleName"]
    df["condition"] = df.mapper.map(patient_map)
    df.condition = df.condition.map(control_condition_rename)
    df["run"] = df.mapper.map(run_map)
    df.rename({"Charge":"charge", "m_score_peptide_global":"searchScore", "Intensity":"intensity", "FullPeptideName":"peptide", "ProteinName":"proteins"},axis=1, inplace = True)
    df.columns
    res = df[["run", "condition", "charge", "searchScore", "intensity", "peptide", "proteins"]]
    return res

def main(filename, mapper_file, output):
    df = convert_pxd021197_to_triqler(filename, mapper_file, output)
    df.to_csv(output, sep = "\t", index = False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script takes triqler results, top3 results, msstat results and msqrob2 results and plots differential HeLa vs differential non-HeLa lineplot.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--filename', type=str,
                        help='PXD021197 input file.')

    parser.add_argument('--mapper_file', type=str,
                        help='PXD021197 annotation file.')

    parser.add_argument('--output', type=str,
                        help='Output file name.')


    # parse arguments from command line
    args = parser.parse_args()
    filename = args.filename
    mapper_file = args.mapper_file
    output = args.output

    main(filename, mapper_file, output)
    
    
#for i in df.condition.unique():
#    print(i, ":", end = " ")
#    print(len(df[df.condition == i].run.unique()))








