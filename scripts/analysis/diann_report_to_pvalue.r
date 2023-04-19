
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library("limma")
library("data.table")
library(tidyverse)
 
design.mat <- cbind(c(1,1,1,0,0,0), c(0,0,0,1,1,1))
contrast.mat <- matrix(c(1,-1), ncol = 1)
dimnames(contrast.mat) <- list(c("A", "B"), "Diff")

sample <- factor(rep(c("A", "B"), each = 3))
design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)

contrast.mat <- makeContrasts(
  Diff = A - B, 
  levels = design.mat)

#df <- fread("report.pg_matrix.tsv", select = c(1,6,7,8,9,10,11))
df <- fread("report.pg_matrix.tsv", select = c(1,6,8,10,7,9,11))
df <- log2(df)
names(df)
df <- df %>% column_to_rownames(., var = 'Protein.Group')


# Fit the expression matrix to a linear model
fit <- lmFit(df, design.mat)
# Compute contrast
fit_contrast <- contrasts.fit(fit, contrast.mat)# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast <- eBayes(fit_contrast)# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast)

head(fit_contrast)

help(topTable)

write.table(fit_contrast$p.value, "limma_p_values.csv", sep = "\t")
write.table(topTable(fit_contrast, n =Inf), "limma_results.csv", sep = "\t")




