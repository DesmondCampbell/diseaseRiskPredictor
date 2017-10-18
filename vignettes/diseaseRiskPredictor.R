## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(diseaseRiskPredictor))

# copy depression disease model into current environment
lDis <- diseaseRiskPredictor::lDisDepression

# report the depression disease model
lDis

## ------------------------------------------------------------------------
# calculate disease model
lDis <- fnCalcDiseaseModel( lDis )
names(lDis)

## ------------------------------------------------------------------------
# copy 3 generation pedigree into current environment
dfPed <- diseaseRiskPredictor::dfPed3Gen1Aff
dfPed

## ------------------------------------------------------------------------
# validate pedigree (also performs type conversion on some columns, hence the returned object)
dfPed <- fnValidatePedigree( dfPed )

## ---- fig.width=7, fig.height=4------------------------------------------
# construct pedigree object
oPedigree <- fnConstructPedigree( dfPed )

# plot pedigree diagram
plot(oPedigree)

## ------------------------------------------------------------------------
lPedDis <- fnPrepareForRiskPrediction( dfPed, lDis )

# calculate pedigree members' risks of depression
lRes <- fnPredictRisk( lPedDis )

# report contents of the results object
names(lRes)

## ---- fig.width=7, fig.height=4------------------------------------------
# report pedigree members' risks
lRes$dfPedRisk[,c(2,3,5:7,11:12)]

# plot pedigree members' risks
fnBarPlotRisk( lRes$dfPedRisk )

## ------------------------------------------------------------------------
lRes$lDis$categoricalRiskFactors[,1:7]

## ---- fig.width=7, fig.height=4------------------------------------------
# report pedigree members' risks
lRes$dfPedRisk[,c(2,3,5:7,13,18)]

## ----eval=FALSE----------------------------------------------------------
#  # install it
#  install.packages("Rcpp")
#  
#  # attach it
#  library(Rcpp)
#  
#  # test c++ code for adding 1 plus 1
#  evalCpp('1+1',showOutput=1,rebuild=1)

## ----eval=FALSE----------------------------------------------------------
#  # install it
#  install.packages("devtools")
#  
#  # attach it
#  library(devtools)
#  
#  # test for Rtools installation
#  find_rtools()

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  
#  devtools::install_github("DesmondCampbell/diseaseRiskPredictor")

## ----eval=FALSE----------------------------------------------------------
#  # attach it
#  library("diseaseRiskPredictor")
#  
#  assignInNamespace( "plot.pedigree", diseaseRiskPredictor:::plot.pedigree.FIXED, ns="kinship2")

