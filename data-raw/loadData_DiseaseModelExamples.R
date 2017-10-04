library(diseaseRiskPredictor)

options(stringsAsFactors = FALSE) # hurray!!!



# change dir to data-raw
sDir <- "data-raw"
if( !grepl( paste0(sDir,"$"),getwd() ) ) setwd(sDir)



sFileDis <- "diseaseModels/dm.depression.txt"
fnRep( 'read disease model from file -', sFileDis )
lDisDepression <- NA; lDisDepression <- fnReadDisFromFile( sFileDis )
fnStr(lDisDepression)

fnRep('save data')
devtools::use_data(lDisDepression, overwrite=T)


sFileDis <- "diseaseModels/dm.schizophrenia.txt"
fnRep( 'read disease model from file -', sFileDis )
lDisSchizophrenia <- NA; lDisSchizophrenia <- fnReadDisFromFile( sFileDis )
fnStr(lDisSchizophrenia)

fnRep('save data')
devtools::use_data(lDisSchizophrenia, overwrite=T)
