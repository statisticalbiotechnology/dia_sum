#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SWATH2stats))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MSstats))
suppressPackageStartupMessages(library("argparse"))

#setwd("/hdd_14T/data/PXD002952/20210805_osw_run")
#filename <- 'concatenated_osw_results_transitions_filtered_n_6.csv'
#output <- "mscore_threshold.csv"
#FFT <- 1
#fdr_target = 0.01

compute_mscore4pepfdr <- function(filename, output, FFT, fdr_target){
  data <- data.frame(fread(filename, sep='\t', header=TRUE))
  Study_design <- data.frame(Filename = unique(data$filename))
  Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
  Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
  Study_design$Run <- seq_len(nrow(Study_design))
  data.annotated <- sample_annotation(data, Study_design, column_file = "filename")
  peptide_filter_m_score_treshold <- mscore4pepfdr(data.annotated, FFT = FFT, fdr_target = fdr_target, mscore.col = "m_score")
  write.table(peptide_filter_m_score_treshold, file = output, sep = "\t", row.names = FALSE, col.names = FALSE) 
  print(paste0("Output: ", output))
}

parser <- ArgumentParser()

parser$add_argument("-i", "--input", 
                    help="Input OSW file. It should be in concatednated format.")
parser$add_argument("-o", "--output", default="mscore_threshold.csv", 
                    help="Print little output")
parser$add_argument("--fft", type="integer", default=1, 
                    help="Ratio of false positives to true negatives, q-values from [Injection_name]_full_stat.csv in pyProphet stats output. As an approximation, the q-values of multiple runs are averaged and supplied as argument FFT. Numeric from 0 to 1. Defaults to 1, the most conservative value (1 Decoy indicates 1 False target). ",
                    metavar="number")
parser$add_argument("--fdr_target", type="double", default=0.01, 
                    help = "FDR target, numeric, defaults to 0.01. An m_score cutoff achieving an FDR < fdr_target will be selected. Calculated as FDR = (TN*FFT/T); TN=decoys, T=targets, FFT=see above. ")

args <- parser$parse_args()

input <- args$input
output <- args$output
fft <- args$fft
fdr_target <- args$fdr_target

compute_mscore4pepfdr(input, output, fft, fdr_target)