library(diseaseRiskPredictor)

options(stringsAsFactors = FALSE) # hurray!!!



# change dir to data-raw
sDir <- "data-raw"
if( !grepl( paste0(sDir,"$"),getwd() ) ) setwd(sDir)


sFilePed <- "pedigrees/ped.3gen.1aff.tsv"
fnRep( paste("Read pedigree from file -", sFilePed ))
# read data from file
dfPed <- NA; dfPed <- read.table( sFilePed, header=T)
fnRep( "Pedigree is" )
dfPed3Gen1Aff <- dfPed

fnRep('save data')
devtools::use_data(dfPed3Gen1Aff, overwrite=T)
