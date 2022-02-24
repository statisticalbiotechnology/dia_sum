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

#data <- data.frame(fread('feature_alignment_filtered_transition.tsv', sep='\t', header=TRUE))
#data <- data.frame(fread('concatenated_osw_output_m_score_0.01.csv', sep='\t', header=TRUE))
#data <- data.frame(fread('concatenated_osw_results.csv', sep='\t', header=TRUE))
#data <- data.frame(fread('osw_output.HYE124_TTOF6600_32fix_lgillet_I150211_002-Pedro_-_Sample_1_-_SW32_-_Repl1.mzML_with_dscore.csv', sep='\t', header=TRUE))
data <- data.frame(fread('concatenated_osw_results_transitions_filtered_n_6.csv', sep='\t', header=TRUE))

Study_design <- data.frame(Filename = unique(data$filename))
Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
Study_design$Run <- seq_len(nrow(Study_design))
head(Study_design)
#write.csv(Study_design, "msstats_run.csv", row.names = FALSE)

data.annotated <- sample_annotation(data, Study_design, column_file = "filename")
#data.annotated.nodecoy <- subset(data.annotated, decoy==FALSE)
rm(data)
gc()

#fdr_target_decoy <- assess_fdr_overall(data.annotated, n.range = 10, 
#                                       FFT = 0.25, output = 'Rconsole')
#mscore4protfdr(data, FFT = 0.25, fdr_target = 0.05)
peptide_filter_m_score_treshold <- mscore4pepfdr(data.annotated, FFT = 1, fdr_target = 0.01, mscore.col = "m_score")

#data.filtered <- filter_mscore_condition(data.annotated, 0.01, n_replica = 2) #0.001 treshold in tutorial
data.filtered <- filter_mscore_condition(data.annotated, mscore = peptide_filter_m_score_treshold, n_replica = 2) #0.001 treshold in tutoria
data.filtered2 <- filter_on_max_peptides(data.filtered, n_peptides = 10)
data.filtered3 <- filter_on_min_peptides(data.filtered2, n_peptides = 2)
rm(data.filtered)
rm(data.filtered2)
rm(data.annotated)
rm(data.annotated.nodecoy)
gc()

saveRDS(data.filtered3, file = "proc_concatenated_osw_results_transitions_filtered_n_6.csv.RData")
data.filtered3 <- readRDS("proc_concatenated_osw_results_transitions_filtered_n_6.csv.RData")
write.csv(data.filtered3, "data.filtered3.csv")
data.transition <- disaggregate(data.filtered3)
MSstats.input <- convert4MSstats(data.transition)

QuantData <- dataProcess(MSstats.input)
comparison <- matrix(c(-1,1), nrow=1)
row.names(comparison) <- "T2-T1"
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)

setwd("/hdd_14T/data/PXD002952/20210805_osw_run/msstats_results")
write.csv(testResultOneComparison$ComparisonResult, "msstat_output_tresholded_20211104.csv", row.names=FALSE)






