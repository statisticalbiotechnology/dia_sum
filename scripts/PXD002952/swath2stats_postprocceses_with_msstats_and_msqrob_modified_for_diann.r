

setwd("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SWATH2stats")


library(SWATH2stats)
library(data.table)
data('Spyogenes', package = 'SWATH2stats')

data_test <- data
data <- data.frame(fread('aligned.csv', sep='\t', header=TRUE))

Study_design <- data.frame(Filename = unique(data$align_origfilename))
Study_design$Condition <- gsub("_.*", "" ,gsub(".*(Sample_)", "", Study_design$Filename))
Study_design$BioReplicate <- gsub("\\..*", "", gsub(".*(Repl)", "", Study_design$Filename))
Study_design$Run <- seq_len(nrow(Study_design))
head(Study_design)

#data.annotated <- sample_annotation(data, Study_design)
data.annotated <- sample_annotation(data, Study_design, column.file = "align_origfilename")
data.annotated.nodecoy <- subset(data.annotated, decoy==FALSE)

count_analytes(data.annotated.nodecoy)

correlation <- plot_correlation_between_samples(data.annotated.nodecoy, column.values = 'Intensity')

correlation <- plot_correlation_between_samples(data.annotated.nodecoy, column.values = 'delta_rt')

variation <- plot_variation(data.annotated.nodecoy)
variation[[2]]

variation_total <- plot_variation_vs_total(data.annotated.nodecoy)
variation_total[[2]]


peptide_signal <- write_matrix_peptides(data.annotated.nodecoy)
protein_signal <- write_matrix_proteins(data.annotated.nodecoy)
head(protein_signal)

par(mfrow = c(1, 3))
fdr_target_decoy <- assess_fdr_overall(data.annotated, n.range = 10, 
                                       FFT = 0.25, output = 'Rconsole')

mscore4protfdr(data, FFT = 0.25, fdr_target = 0.05)

data.filtered <- filter_mscore_condition(data.annotated, 0.001, n.replica = 2)
data.filtered2 <- filter_on_max_peptides(data.filtered, n_peptides = 10)
data.filtered3 <- filter_on_min_peptides(data.filtered2, n_peptides = 2)
data.transition <- disaggregate(data.filtered3)

MSstats.input <- convert4MSstats(data.transition)
MSstats.input <- convert4MSstats(disaggregate(data.annotated)) # NOTE HTIS MODIFICATION TO UNFILTERED DATA
head(MSstats.input)
write.csv(MSstats.input, "msstats_20210511_unfiltered.csv", row.names = FALSE)

mapDIA.input <- convert4mapDIA(data.transition)
head(mapDIA.input)
write.csv(mapDIA.input, "mapDIA.csv", row.names = FALSE)

aLFQ.input <- convert4aLFQ(data.transition)
head(aLFQ.input)

sessionInfo()




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSstats")
library(MSstats)



AADGSTVAQTALSYDDYR_2_44730_y11_1_NA
AADGSTVAQTALSYDDYR

MSstats.input[MSstats.input$PeptideSequence == "AADGSTVAQTALSYDDYR",]
MSstats.input[MSstats.input$PeptideSequence == "VDVGQQPLR",]
QuantData[QuantData$FEATURE == "AADGSTVAQTALSYDDYR_2_44730_y11_1_NA", ]
qMat <- QuantData[1]
qMat <- data.frame(qMat)
qMat[qMat$ProcessedData.FEATURE == "AADGSTVAQTALSYDDYR_2_44730_y11_1_NA", ]
QuantData <- dataProcess(MSstats.input)

write.csv(QuantData, "msstat_output_20210512.csv", row.names = FALSE)

#### TEST with OSW ######
setwd("/hdd_14T/data/PXD002952/osw_res_20210303/hye124/ttof6600/32fix")

MSstats.input <- read.csv("msstats_20210511_unfiltered.csv")



#### START HERE ##### 20210603
setwd("/hdd_14T/data/PXD002952/res_20210530_DIAUmpire/MSFragger")

MSstats.input.diann <- read.csv("msstats_input_diann.tsv")

QuantData <- dataProcess(MSstats.input)

#df[df$aged <= df$laclen, ] 
QuantData
dataProcessPlots(data=QuantData, type="ProfilePlot")

# Quality control plot
dataProcessPlots(data=QuantData, type="QCPlot")	

#  # Quantification plot for conditions
dataProcessPlots(data=QuantData, type="ConditionPlot")

levels(QuantData$ProcessedData$GROUP_ORIGINAL)
comparison <- matrix(c(-1,1), nrow=1)
row.names(comparison) <- "T2-T1"

# Tests for differentially abundant proteins with models:
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)
testResultOneComparison$ComparisonResult
typeof(testResultOneComparison$ComparisonResult)
testResultOneComparison$ComparisonResult

#####################
##### MSqRobSum #####
#####################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MSnbase")

library(MSnbase)
install.packages("reticulate")
install.packages("ncdf4")



## Location of the Maxquant output file
data_path = system.file('extdata','peptides.txt.gz', package = 'msqrobsum')
## file is compressed with gzip (hence '.gz' extention), thus first make connection
data_path = gzfile(data_path)
## Look for the columns with peptide intensities in the maxquant output (columns starting with `Intensity `, 1 column for every sample.)
exprs_col = grepEcols(data_path, 'Intensity ',split = '\t')
## read data into a MSnSet object. The different peptides are indicated in the `Sequence` column.
## All other columns were added as a featureData dataframe in the MSnSet object 
set = readMSnSet2(data_path ,ecol = exprs_col,fnames = 'Sequence'
                  , sep = '\t',stringsAsFactors = FALSE)


library(MSqRobSum)

results1 <- msqrobsum( data = peptide_intensities
                       , mode = 'sum'
                       , group_vars = 'protein'
                       , parallel_args = list(strategy = 'sequential'))


