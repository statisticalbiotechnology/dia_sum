
suppressPackageStartupMessages(library(SWATH2stats))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MSstats))
suppressPackageStartupMessages(library("argparse"))


input <- "Ctrl_ST_msstats_input.csv"
output <- "Ctrl_ST_msstats_results.csv"
output_protein <- "CTRL_ST_msstats_protein_results.csv"

input <- "Ctrl_LT_msstats_input.csv"
output <- "Ctrl_LT_msstats_results.csv"
output_protein <- "CTRL_LT_msstats_protein_results.csv"

input <- "LT_ST_msstats_input.csv"
output <- "LT_ST_msstats_results.csv"
output_protein <- "LT_ST_msstats_protein_results.csv"

run_msstats <- function(input, output, output_protein){
  data <- data.frame(fread(input, sep = ',', header = TRUE))
  MSstats.input <- data
  QuantData <- dataProcess(MSstats.input)
  write.csv(QuantData$ProteinLevelData, output_protein, row.names=FALSE)
  comparison <- matrix(c(-1,1), nrow=1)
  row.names(comparison) <- "T2-T1"
  colnames(comparison) <- c("ST","LT") #MODIFY HERE
  testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData)
  write.csv(testResultOneComparison$ComparisonResult, output, row.names=FALSE)
  
}

parser <- ArgumentParser()

parser$add_argument("--input", 
                    help = "MSstats input file. NOTE: The input file needs to be formatted and filtered before using this script.")
parser$add_argument("--output", 
                    help = "Output name.")
parser$add_argument("--output_protein", 
                    help = "Protein quantity output file name.")
args <- parser$parse_args()

run_msstats(input = args$input, output = args$output, output_protein = args$output_protein)
