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
#  devtools::install_github("diseaseRiskPredictor")
#  
#  # attach it
#  library("diseaseRiskPredictor")

## ----eval=FALSE----------------------------------------------------------
#  assignInNamespace( "plot.pedigree", diseaseRiskPredictor:::plot.pedigree.FIXED, ns="kinship2")

