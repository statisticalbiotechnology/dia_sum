
library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(gridExtra)


# file firectory 
#os.chdir("/hdd_14T/data/PXD002952/20210805_osw_run")
setwd("/hdd_14T/data/PXD002952/20210805_osw_run")
peptidesFile = ("/hdd_14T/data/PXD002952/20210805_osw_run/msqrob2_input_20220204.tsv")



ecols <- grep("HYE124\\_", names(read.delim(peptidesFile)))

pe <- readQFeatures(
  table = peptidesFile, fnames = 1, ecol = ecols,
  name = "peptideRaw", sep = "\t"
)

cond <- 1

colData(pe)$condition <- substr(colnames(pe), cond, cond) %>%
  unlist() %>%
  as.factor()


rowData(pe[["peptideRaw"]])$nNonZero <- rowSums(assay(pe[["peptideRaw"]]) > 0)

pe <- zeroIsNA(pe, "peptideRaw") # convert 0 to NA



MSnbase::plotNA(assay(pe[["peptideRaw"]])) +
  xlab("Peptide index (ordered by data completeness)")

pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
limma::plotDensities(assay(pe[["peptideLog"]]))

pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins))

pe <- filterFeatures(pe, ~ nNonZero >= 2)
nrow(pe[["peptideLog"]])

pe <- normalize(pe,
                i = "peptideLog",
                name = "peptideNorm",
                method = "center.median"
)

limma::plotDensities(assay(pe[["peptideNorm"]]))


boxplot(assay(pe[["peptideNorm"]]),
        col = palette()[-1],
        main = "Peptide distribtutions after normalisation", ylab = "intensity"
)

limma::plotMDS(assay(pe[["peptideNorm"]]), col = as.numeric(colData(pe)$condition))

pe <- aggregateFeatures(pe,
                        i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                        name = "protein"
)


plotMDS(assay(pe[["protein"]]), col = as.numeric(colData(pe)$condition))


pe <- msqrob(object = pe, i = "protein", formula = ~condition)

getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])

L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
pe <- hypothesisTest(object = pe, i = "protein", contrast = L)


msqrob2_df <-rowData(pe[["protein"]])$conditionB

write.csv(msqrob2_df, "msqrob2_results.tsv")
