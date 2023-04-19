setwd("/hdd_14T/data/PXD002952/20210805_osw_run")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SWATH2stats")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")

library(SWATH2stats)
library(data.table)
library(MSstats)

data <- data.frame(fread('concatenated_osw_results.csv', sep='\t', header=TRUE))


Study_design <- data.frame(Filename = unique(data$filename))
Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
Study_design$Run <- seq_len(nrow(Study_design))
head(Study_design)

data.annotated <- sample_annotation(data, Study_design, column_file = "filename")
rm(data)

#mscore4protfdr(data, FFT = 0.25, fdr_target = 0.05)
peptide_filter_m_score_treshold <- mscore4pepfdr(data.annotated, FFT = 1, fdr_target = 0.01, mscore.col = "m_score")
data.filtered <- filter_mscore_condition(data.annotated, mscore = peptide_filter_m_score_treshold, n_replica = 2) #0.001 treshold in tutoria


write.csv(data.filtered,"peptide_filtered_osw_result.csv", row.names = FALSE)





