setwd("/hdd_14T/data/PXD002952/20210614_dataset/diaumpire_spectral_lib_20210706/MSFragger_20210707/diann_20210811")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SWATH2stats")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MSstats")

library(SWATH2stats)
library(data.table)
library(MSstats)

data <- data.frame(fread('diann_msstats_input.csv', sep='\t', header=TRUE))
#data <- data.frame(fread('diann_msstats_input_recomputed_fdr_20211201.csv', sep = ',', header = TRUE))

Study_design <- data.frame(Filename = unique(data$filename))
Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
Study_design$Run <- seq_len(nrow(Study_design))
head(Study_design)
write.csv(Study_design, "msstats_run.csv", row.names = FALSE)

data.annotated <- sample_annotation(data, Study_design, column.file = "align_origfilename")
data.annotated.nodecoy <- subset(data.annotated, decoy==FALSE)

fdr_target_decoy <- assess_fdr_overall(data.annotated, n.range = 10, 
                                       FFT = 0.25, output = 'Rconsole')
mscore4protfdr(data, FFT = 0.25, fdr_target = 0.05)
data.filtered <- filter_mscore_condition(data.annotated, 0.001, n.replica = 2) #0.001 treshold in tutorial
data.filtered2 <- filter_on_max_peptides(data.filtered, n_peptides = 10)
data.filtered3 <- filter_on_min_peptides(data.filtered2, n_peptides = 2)
data.transition <- disaggregate(data.filtered3)
MSstats.input <- convert4MSstats(data.transition)

MSstats.input <- data
QuantData <- dataProcess(MSstats.input)
comparison <- matrix(c(-1,1), nrow=1)
row.names(comparison) <- "T2-T1"
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)

#setwd("/hdd_14T/data/PXD002952/20210805_osw_run/msstats_results")
write.csv(testResultOneComparison$ComparisonResult, "msstat_output_tresholded_20211104.csv", row.names=FALSE)
