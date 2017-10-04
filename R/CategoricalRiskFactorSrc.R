# This code implements a solution to the following problem
#
# Problem Statement
# We wish to convert epi info on a categorical risk factor for a dichotomous trait into liability distributions per risk factor stratum.
# The underlying principle we will use to do that is the liability threshold model.
# We have the following info
#   Freqs of each risk factor stratum
#   The relative risk of each stratum relative to a reference stratum
#   trait lifetime risk
# We will assume disease liability within a stratum is Gaussian and the within stratum liability variance is the same across all strata.
#
# D Campbell 26.5.15



#' Validates categorical risk factor inputs, calculates directly computable parameters
#'
#' Validates categorical risk factor inputs are in correct ranges etc.
#' Calculates directly computable risk factor parameters.
#' Checks whether these are valid - some input combos are invalid, e.g. produce stratum risks > 1.
#' @param lifetimeRisk population lifetime risk for the disease
#' @param vFreq population frequency vector for strata (must sum to 1)
#' @param vRelRisk relative risk vector for strata (one element must equal 1, this identifies the reference stratum)
#' @return data.frame containing directly calculable risk factor parameters
#' @export
fnValidateCategoricalRiskFactorInputs <- function( lifetimeRisk, vFreq, vRelRisk )
{
  # validate inputs

  # validate lifetime risk
  if( lifetimeRisk<=0 | lifetimeRisk>=1 ) stop("Disease lifetime risk must be between 0 and 1") # validate lifetime risk

  # validate risk strata freqs
  if( any(vFreq<=0) | any(vFreq>=1) ) stop("risk strata frequencies must be between 0 and 1")

  # validate relative risks
  if( any(vRelRisk<=0 )) stop("Relative risks must be greater than 0")

  if(length(vFreq) != length(vRelRisk)) stop("Expected vectors vFreq and vRelRisk to have the same length")

  # validate more than one stratum
  if(length(vFreq) <= 1 ) stop("Categorical risk factor must have more than one stratum")

  # init obj
  dfRF <- data.frame( freq=vFreq, relRisk=vRelRisk)
  dfRF

  ## validate inputs combo

  # validate risk strata freqs sum to 1
  if( abs(sum(dfRF$freq)-1)>1e-4 ) stop("risk strata frequencies must be sum to 1, but in fact sum to ", sum(dfRF$freq) )
  # if sum differs from 1 by rounding error adjust so as to sum to 1
  dfRF$freq <- dfRF$freq / sum(dfRF$freq)

  # identify a ref stratum
  ixRefStratum <- which(dfRF$relRisk == 1)[1]
  if(is.na(ixRefStratum)) stop("One risk factor stratum must be specified as reference stratum, by giving it a relative risk of 1")
  dfRF$bRefStratum <- F; dfRF$bRefStratum[ixRefStratum] <- T

  vRelRiskNonRef <- dfRF$relRisk[-ixRefStratum]
  vStrataFreqNonRef <- dfRF$freq[-ixRefStratum]

  # calc reference stratum freq and risk
  f0 <- 1 - sum(vStrataFreqNonRef)
  r0 <- lifetimeRisk/(f0+sum(vRelRiskNonRef*vStrataFreqNonRef))

  # calc risks per stratum
  dfRF$risk <- r0*dfRF$relRisk

  # validate risks
  if( any(dfRF$risk>=1 | any(dfRF$risk<=0) ) ){
    sErr <- "The categorical risk factor info and disease lifetime risk imply impossible per stratum lifetime risks (see above)"
    print(dfRF)
    stop(sErr)
  }

  # calc proportion of population affected
  dfRF$propPopAff <- dfRF$freq*dfRF$risk

  list( lifetimeRisk=lifetimeRisk, dfRF=dfRF )
}



fnOptimDiscrepancy <- function( vLiaMeanNonRef, vStrataFreqNonRef, lifetimeRisk, vRelRiskNonRef )
{
  #### validate inputs
  if( any(vRelRiskNonRef<=0)) stop("relative risks must be positive non-zero finite")

  if( lifetimeRisk<0) stop("lifetimeRisk invalid (<0)")
  if( lifetimeRisk>1) stop("lifetimeRisk invalid (>1)")

  if( any(vStrataFreqNonRef<0)) stop("Non Ref Stratum freqs contains 1 or more bad freqs (<0)")
  if( any(vStrataFreqNonRef>1)) stop("Non Ref Stratum freqs contains 1 or more bad freqs (>1)")
  if( sum(vStrataFreqNonRef)>=1) stop("Non Ref Stratum freqs must sum to less than 1")


  #### calc parameters whose value is implied by the known fixed inputs

  # calc ref stratum freq
  f0 <- 1 - sum(vStrataFreqNonRef)
  if( f0 == 0 ) stop("Reference Stratum freq must be non-zero")

  # calc ref stratum risk
  r0 <- lifetimeRisk/(f0+sum(vRelRiskNonRef*vStrataFreqNonRef))

  # calc expected lifetime risks
  # note inside this function ref stratum is always the first stratum, that is not necessarily the case outside of it
  vLRiskExp <- r0*c(1,vRelRiskNonRef)
  #### validate
  if( any(vLRiskExp>=1) ) stop("The categorical risk factor info and disease lifetime risk imply impossible per stratum lifetime risks")

  # calc ref stratum liability mean
  u0 <- -sum(vLiaMeanNonRef*vStrataFreqNonRef)/f0


  vLiaMean <- c(u0,vLiaMeanNonRef)
  vStrataFreq <- c(f0,vStrataFreqNonRef)

  # calc liability variance explained by the risk factor
  liaVarInter <- sum(vStrataFreq*(vLiaMean^2))

  # calc liability variance not explained by the risk factor
  liaVarIntra <- 1 - liaVarInter

  # return a large discrepancy when the current solution implies the intra statum variance is negative
  if(liaVarIntra<0) {
    # the max possible dev of risk is 1
    # riskDev = sum(vStrataFreq*(vLRiskImp - vLRiskExp)^2) = sum(vStrataFreq*(maxDevRisk)^2) = sum(vStrataFreq*(1)^2) = 1
    # increment this with the negative variance to give it directionality for algorithms that require gradient
    riskDev <- 1-liaVarIntra
    dfLog <<- rbind(dfLog, data.frame(u0=u0, liaVarIntra=liaVarIntra, liaThr=NA, liaMeanImp=NA, riskDev=riskDev ))
    return(riskDev)
  }

  # calc liability threshold for ref stratum
  liaSdIntra <- sqrt(liaVarIntra)
  liaThr <- stats::qnorm(1-r0, mean=u0, sd=liaSdIntra )

  # calc implied lifetime risks
  vLRiskImp <- stats::pnorm( q = liaThr, mean = vLiaMean, sd = liaSdIntra, lower.tail = F)

  # calc discrepancy

  # the following weights vector is appropraite for the average discrepancy in the population
  #vWeights <- vStrataFreq
  # the following weights vector gives equal weighing to the discrepancy in the population and the discrepancy across strata
  # it upweighs very rare strata, this guards against the possibility that a good fit for the population may come at the price of a disastrous fit for the very rare
  vWeights <- (vStrataFreq + 1/length(vStrataFreq))
  vWeights <- vWeights/sum(vWeights)

  # calc discrepancy - weighted distance between the guessed lifetime risks and the actual
  riskDev <- sum(vWeights*(vLRiskImp - vLRiskExp)^2)

  liaMeanImp <- sum(vStrataFreq*vLiaMean)

  dfLog <<- rbind(dfLog, data.frame(u0=u0, liaVarIntra=liaVarIntra, liaThr=liaThr, liaMeanImp=liaMeanImp, riskDev=riskDev ))

  # return the discrepancy
  riskDev
}



#' Optimises liability mean and variance for a categorical risk factor
#'
#' Uses optimisation to find the liability mean and variance for the strata of a categorical risk factor
#' @param lifetimeRisk population lifetime risk for the disease
#' @param dfRF directly calculable risk factor parameters
#' @param sMethod optimisation method
#' @param lControl optimisation method tuning parameters
#' @return list containing solution
#' @export
fnRunOptim <- function( lifetimeRisk, dfRF, sMethod, lControl )
{
  lResult <- list( lifetimeRisk=lifetimeRisk, dfRF=dfRF )

  #### initialise optimiser
  vLiaMeanNonRef <- rep(0,nrow(dfRF)-1) # initial guess
  vStrataFreqNonRef <- dfRF$freq[dfRF$bRef==F]
  vRelRiskNonRef <- dfRF$relRisk[dfRF$bRef==F]
  dfLog <<- NULL # this global must be reset to null before calling optim()
  # optimise
  lOptim <- NA; lOptim <- stats::optim( vLiaMeanNonRef, fnOptimDiscrepancy, vStrataFreqNonRef=vStrataFreqNonRef, lifetimeRisk=lifetimeRisk, vRelRiskNonRef=vRelRiskNonRef, control=lControl, method=sMethod )

  # add optimiser settings and log to results
  lOptim$lControl <- lControl
  lOptim$sMethod <- sMethod
  lOptim$dfLog <- dfLog

  lResult$lOptim <- lOptim

  # check whether optimiser converged on a solution
  if( lOptim$convergence ) {
    lResult$sErrMsg <- "optimiser failed to converge"
    return(lResult)
  }

  # convert solution

  # find ref stratum
  ixRefStratum <- which(dfRF$bRefStratum)
  if(length(ixRefStratum)!=1) stop("Expected exactly 1 ref stratum")

  # calc per strata lia distribution
  dfRF$liaMean <- NA
  dfRF$liaMean[-ixRefStratum] <- lOptim$par
  dfRF$liaMean[ixRefStratum] <- -sum(dfRF$liaMean[-ixRefStratum]*dfRF$freq[-ixRefStratum])/dfRF$freq[ixRefStratum]

  # calc inter stratum variance
  liaVarInter <- sum(dfRF$freq*dfRF$liaMean^2)
  dfRF$liaVarInter <- liaVarInter

  lResult$dfRF <- dfRF # update it

  # calc intra stratum variance
  liaVarIntra <- 1-liaVarInter
  if(liaVarIntra<0) {
    lResult$sErrMsg <- "Solution implies negative intra stratum variance"
    return(lResult)
  }
  dfRF$liaVarIntra <- liaVarIntra

  # calc thresholds implied by solution
  dfRF$liaThr_imp <- stats::qnorm(dfRF$risk,dfRF$liaMean,sqrt(dfRF$liaVarIntra),lower.tail=F)

  # calc risks implied by solution (using ref stratum lia threshold, as this is the threshold the discrepancy function uses)
  liaThr0 <- dfRF$liaThr_imp[ixRefStratum]
  dfRF$risk_imp <- stats::pnorm( liaThr0, dfRF$liaMean, sqrt(dfRF$liaVarIntra), lower.tail=F)

  # calc proportion of population affected implied by solution
  dfRF$propPopAff_imp <- dfRF$freq*dfRF$risk_imp

  lResult$dfRF <- dfRF # update it

  # calc quality of solution
  dfRFDiagnostics <- data.frame(
    liaThr_imp_var = stats::var(dfRF$liaThr_imp),
    risk_imp_mse=mean((dfRF$risk_imp-dfRF$risk)^2),
    popLiaMean_dev=sum(dfRF$liaMean*dfRF$freq),
    lifetimeRisk_dev=sum(dfRF$propPopAff_imp)-lifetimeRisk )
  dfRFDiagnostics

  lResult$dfRFDiagnostics <- dfRFDiagnostics

  # test for a bad solution
  sErrMsg <- NULL
  if( dfRFDiagnostics$liaThr_imp_var > 0.01 ) sErrMsg <- c( sErrMsg, "Variance of the implied liability thresholds over strata exceeds limit" )
  if( dfRFDiagnostics$risk_imp_mse > 0.01 ) sErrMsg <- c( sErrMsg, "Mean Squared Deviation of the implied risks over strata exceeds limit" )
  if( abs(dfRFDiagnostics$popLiaMean_dev) > 0.01 ) sErrMsg <- c( sErrMsg, "Deviation of the implied population mean from zero exceeds limit" )
  if( abs(dfRFDiagnostics$lifetimeRisk_dev) > 0.01 ) sErrMsg <- c( sErrMsg, "Deviation of the implied lifetime risk from true value exceeds limit" )

  lResult$sErrMsg <- sErrMsg
  lResult
}



#' Generates plots indicative of whether Categorical risk factor optimisation converged
#'
#' Generates 2 plots x axis being number of iterations, y axis being
#' - categorical risk factor reference stratum mean
#' - discrepancy
#' The reference stratum mean should stabalise.
#' the discrepancy should fall to a very low value.
#' @param sMethod optimisation method used
#' @param args categorical risk factor info
#' @param dfLog log containing history of optimisation
#' @return NULL
#' @export
plotRefStraMean <- function(sMethod, args, dfLog)
{
  dfLog$iter <- 1:nrow(dfLog)
  iDropped <- 10
  if( nrow(dfLog) > 10*iDropped ) dfLog <- dfLog[-(1:iDropped),]

  sSuffix <- paste( paste(args, collapse=", "), sMethod, sep="\n")
  par(mfrow=c(1,2))
  sMain <- paste("Ref Stratum Lia Mean\n", sSuffix )
  suppressWarnings( plot(x=dfLog$iter, y=dfLog$u0, main=sMain))
  sMain <- paste("Discrepancy that drives optimisation\n", sSuffix )
  suppressWarnings( plot(x=dfLog$iter, y=dfLog$riskDev,log="y", main=sMain))
}
