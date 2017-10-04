library(diseaseRiskPredictor)

options(stringsAsFactors = FALSE) # hurray!!!

# change dir to data-raw
sDir <- "data-raw"
if( !grepl( paste0(sDir,"$"),getwd() ) ) setwd(sDir)

sFileDialect <- "disFileParameters4.tsv"
dfDialectDisFile <- read.table(sFileDialect, sep="\t",header=T)
fnStr(dfDialectDisFile)



fnRep('save data')
devtools::use_data(dfDialectDisFile, overwrite=T)
