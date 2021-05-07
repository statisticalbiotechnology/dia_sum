
#Check msqrobsum.html in github repo for more information...

expr_path <- "expr.csv"
fd_path <- "fd.csv"
pd_path <- "pd.csv"

exprs_col = grepEcols(expr_path, '00',split = '\t')

set = readMSnSet2(expr_path ,ecol = exprs_col,fnames = 'peptide'
                  , sep = '\t',stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(tidyverse))

expr <- read.csv2("expr.csv", sep = "\t", row.names = 1)
fd <- read.csv("fd.csv", sep = "\t", row.names = 1)
pd <- read.csv("pd.csv", sep = "\t", row.names = 1)
fData(set) = fd
pData(set) = pd
exprs(set)

exprs(set)[0 == (exprs(set))] <- NA
set <- set[!fData(set)$decoy]
set_log = log(set, base = 2)
set = normalize(set, 'vsn') #vsn normalized values are on log-scale


library(limma)
plotDensities(exprs(set_log))
plotMDS(exprs(set_log), top = Inf,col = as.integer(pData(set)$condition))
plotDensities(exprs(set))
plotMDS(exprs(set), top = Inf,col = as.integer(pData(set)$condition))

protset <- suppressWarnings(combineFeatures(set,fun="robust", groupBy = fData(set)$protein))

library(msqrobsum)
set_df = MSnSet2df(set)
set_df

protset_df = MSnSet2df(protset)
protset_df

#msprobsum on peptide intentities.
msqrobsum_result <- msqrobsum(data = set, expression ~ (1 | condition)
                              , contrasts = 'condition', mode = 'msqrob')

saveRDS(msqrobsum_result, file = "msrobsum_results.rds")

#msqprobsum on protein summarize.
msqrobsum_result <- msqrobsum(data = protset, expression ~ (1 | condition)
                              , contrasts = 'condition', mode = 'msqrobsum'
                              ## group by folowing variables,
                              ## they will also be retained in output
                              , group_vars = c('protein','human','ecoli'))
msqrobsum_result


# USE THIS!!!

formulas =  c(expression ~ (1|condition) + (1|sample) + (1|feature)
              , expression ~ (1|condition))

msqrob_result <- msqrobsum(data = set, formulas, contrasts = 'condition', mode = 'msqrob'
                           ## group by folowing variables,
                           ## they will also be retained in output
                           , group_vars = c('protein','human','yeas8','ecoli'))

saveRDS(msqrobsum_result, file = "msrobsum_results.rds")

typeof(msqrobsum_result)
class(msqrobsum_result)

msqrobsum_result
  
