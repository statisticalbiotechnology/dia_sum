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

data <- data.frame(fread('osw_msstats_input.csv', sep = ',', header = TRUE))
MSstats.input <- data
QuantData <- dataProcess(MSstats.input)
comparison <- matrix(c(-1,1), nrow=1)
row.names(comparison) <- "T2-T1"
colnames(comparison) <- c(1,2)
testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)

#setwd("/hdd_14T/data/PXD002952/20210805_osw_run/msstats_results")
write.csv(testResultOneComparison$ComparisonResult, "msstat_output_tresholded_20222121.csv", row.names=FALSE)
