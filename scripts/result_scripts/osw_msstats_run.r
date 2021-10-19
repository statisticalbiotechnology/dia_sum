library(SWATH2stats)
library(data.table)


setwd("/hdd_14T/data/PXD002952/20210805_osw_run")

#data <- data.frame(fread('feature_alignment.tsv', sep='\t', header=TRUE))
data <- data.frame(fread('feature_alignment_filtered_transition.tsv', sep='\t', header=TRUE))

Study_design <- data.frame(Filename = unique(data$align_origfilename))
Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
Study_design$Run <- seq_len(nrow(Study_design))

data.annotated <- sample_annotation(data, Study_design, column.file = "align_origfilename")
data.annotated.nodecoy <- subset(data.annotated, decoy==FALSE)

#data.filtered <- filter_mscore_condition(data.annotated, 0.01, n.replica = 2)
#data.filtered2 <- filter_on_max_peptides(data.filtered, n_peptides = 10)
#data.filtered3 <- filter_on_min_peptides(data.filtered2, n_peptides = 2)
#data.transition <- disaggregate(data.filtered3)
#MSstats.input <- convert4MSstats(data.transition)

MSstats.input <- convert4MSstats(disaggregate(data.annotated)) # NOTE HTIS MODIFICATION TO UNFILTERED DATA



