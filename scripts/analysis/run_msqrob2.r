suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(QFeatures))
suppressPackageStartupMessages(library(msqrob2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library("argparse"))

run_msqrob2 <- function(input, output, output_protein){
  peptidesFile = (input)
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
  pe <- logTransform(pe, base = 2, i = "peptideRaw", name = "peptideLog")
  
  pe <- filterFeatures(pe, ~ Proteins %in% smallestUniqueGroups(rowData(pe[["peptideLog"]])$Proteins))
  pe <- filterFeatures(pe, ~ nNonZero >= 2)
  pe <- normalize(pe,
                  i = "peptideLog",
                  name = "peptideNorm",
                  method = "center.median"
  )
  
  pe <- aggregateFeatures(pe,
                          i = "peptideNorm", fcol = "Proteins", na.rm = TRUE,
                          name = "protein"
  )
  
  pe <- msqrob(object = pe, i = "protein", formula = ~condition)
  getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
  L <- makeContrast("conditionB=0", parameterNames = c("conditionB"))
  pe <- hypothesisTest(object = pe, i = "protein", contrast = L)
  
  msqrob2_df <-rowData(pe[["protein"]])$conditionB
  write.csv(msqrob2_df, output)
  write.csv(assay(pe[["protein"]]), output_protein)
}


parser <- ArgumentParser()

parser$add_argument("--input", 
                    help = "MSqRob2 input file. NOTE: The input file needs to be formatted and filtered before using this script.")
parser$add_argument("--output", 
                    help = "Output name.")
parser$add_argument("--output_protein", 
                    help = "Output name for protein quantity file.")
args <- parser$parse_args()

run_msqrob2(input = args$input, output = args$output, output_protein = args$output_protein)




