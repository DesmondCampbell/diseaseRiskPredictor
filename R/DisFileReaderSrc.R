# This is the cost function for an optimisation routine.
# The optimisation routine calls this function repeatedly adjusting the function inputs each time, so as to minimise the value this function returns.
# The returned value is the discrepancy between the true lifetime risk vector and the lifetime risk vector implied by the supplied liability mean vector.
# The discrepancy is quantified as a Euclidean distance.
#
# XXXX BUG - I have found the optim function in combination with this error function yields solutions which closely match the desired relative risks.
# The solution should also be such that the mean liability aggregating across all strata is zero. This constraint is not imposed, despite this the solutions found do pretty well match this constraint.
# However it would be even better if we imposed the constraint, i.e. vStrataFreq*vLiaMean=0.
# Ways to do this suggest themselves
# 1) impose the vStrataFreq*vLiaMean=0 via the constOptim function
# 2) add the following to the value returned by this function
#         (sum(vStrataFreq*vLiaMean))^2
# 3) instead of supplying optim() and fnCalcError() with vLiaMean, supply them with vLiaMean[-1] then calculate the missing vLiaMean[1] under the constraint
#
fnCalcError <- function( vLiaMean, vStrataFreq, liaThr, vLRisk  )
{
  # calc liability variance explained by the risk factor
  liaVarInter <- sum(vStrataFreq*(vLiaMean^2))
  # calc liability variance not explained by the risk factor
  liaVarIntra <- 1 - liaVarInter

  # XXXX
  # given a bad guess of vLiaMean, liaVarIntra can be negative which means liaSdIntra can't be calculated, in which circumstances the program breaks
  # normally would use constrOptim() to constrain vLiaMean to always be sensible, but it can't be used here because the appropriate constraint can't be specified
  # so instead when we have a negative variance, we set std dev to an appropriate value
  liaSdIntra <- NA; if(liaVarIntra>0) {
    liaSdIntra <- sqrt(liaVarIntra)
  } else {
    liaSdIntra <- 10^(liaVarIntra-2)
  }

  # calc implied lifetime risk
  vLRiskGuess <- stats::pnorm( q=liaThr, mean = vLiaMean, sd = liaSdIntra, lower.tail = F)

  # calc error - the difference between the guessed lifetime risks and the actual
  sum(vStrataFreq*(vLRiskGuess - vLRisk)^2)
}



fnFindRiskFactors <- function( lDis )
{
  # find all risk factors
  vsuCatRF <- NULL
  if( "categoricalRiskFactors" %in% names(lDis)) {
    vsuCatRF <- unique(lDis$categoricalRiskFactors$riskFactor)
  }
  vsuCatRF
  vsuCovRF <- NULL
  vsPartition <- c("Asq","Csq","Esq")
  for( sPart in c("Tot",vsPartition) ){
    sMx <- paste0("covMatrix",sPart)
    if(!is.null(lDis[[sMx]])){

      vsuCovRF <- c( vsuCovRF, rownames(lDis[[sMx]]))
      vsuCovRF <- c( vsuCovRF, colnames(lDis[[sMx]]))
      vsuCovRF <- unique( vsuCovRF )
    }
  }
  vsuAllRF <- unique(c("liability",vsuCatRF,vsuCovRF))

  # find quantitative risk factors
  vsuQuaRF <- setdiff(vsuCovRF, c("liability",vsuCatRF))
  vsuQuaRF
  vsuCatRF

  lDis$vsuAllRF <- vsuAllRF
  lDis$vsuCatRF <- vsuCatRF
  lDis$vsuQuaRF <- vsuQuaRF
  lDis
}



#' Reads a disease model file into a named list of values
#' XXXX
#' @param sFileDis disease model file
#' @return lDis disease model
#' @export
fnReadDisFromFile <- function( sFileDis ){

  # read file into a vector of lines
  vsLines <- readLines( sFileDis, n = -1)

  lDis <- fnReadLinesIntoDisList(vsLines)
  lDis
}


#' Reads a disease model from file
#'
#' Parses the supplied disease model file according to file format defined in the supplied disease model file format definition.
#' @param sFileDis path to disease model file
#' @return disease model
#' @export
fnReadFileDis <- function( sFileDis )
{
  fnRep("fnReadFileDis()")

  # read file into a vector of lines
  vsLines <- readLines( sFileDis, n = -1)

  # construct lDis object from file contents
  lDis <- fnReadLinesDis( vsLines )

  lDis
}



#' XXXX
#' XXXX
#' @param vsLines XXXX
#' @return lDis XXXX
#' @export
fnReadLinesIntoDisList <- function( vsLines )
{
  # read definition for Dis file
  dfDisDef <- diseaseRiskPredictor::dfDialectDisFile

  # extract all possible Dis file parameters from Dis file as text
  vbx <- dfDisDef$category == "global" & dfDisDef$type %in% c("character","characterVector")
  vsTokenLinesAfter <- dfDisDef$parameter[vbx]

  vbx <- dfDisDef$category == "global" & ! dfDisDef$type %in% c("character","characterVector","table")
  vsTokenNameValuePair <- dfDisDef$parameter[vbx] # XXXX what about SUFFIX=beta
  # XXXX vsTokenNameValuePair <- c("heritability","sharedEnvironmentability","lifetimeRisk")

  vbx <- dfDisDef$category == "global" & dfDisDef$type %in% c("table")
  vsTokenTableAfter <- dfDisDef$parameter[vbx]

  # convert into read object
  lDis <- fnReadLines( vsLines, vsTokenNameValuePair, vsTokenLinesAfter, vsTokenTableAfter )

  # convert parameter values from text to appropriate type and do parameter value specific validation
  lDis <- fnTypeConvertRangeCheck( dfDisDef, lDis )

  lDis
}



#' Calculates disease model based on epidemiological findings
#'
#' Calculates liability threshold based disease model based on epidemiological findings
#' @param lDis list containing disease epidemiological findings
#' @return lDis list containing disease epidemiological findings and disease model parameters
#' @export
fnCalcDiseaseModel <- function( lDis )
{
  # read definition for Dis file
  dfDisDef <- diseaseRiskPredictor::dfDialectDisFile

  vsProblem <- NULL
  lDisTmp <- tryCatch( {

    # do extra validation specific to dis

    # add ons
    lDis$liaThr <- stats::qnorm( p=lDis$lifetimeRisk, mean = 0, sd = 1, lower.tail = F)

    # convert the data frames read from the file into matrices
    # label the rows of the cov matrix with the id column values then drop the id column from the matrix
    vsPartition <- c("Asq","Csq","Esq")
    for( sPart in c("Tot",vsPartition) ){
      sMx <- paste0("covMatrix",sPart)
      if(!is.null(lDis[[sMx]])){

        rownames(lDis[[sMx]]) <- lDis[[sMx]]$id
        lDis[[sMx]]$id <- NULL

        # replace by matrix of the same
        vX <- apply( lDis[[sMx]], 2, as.numeric )
        mMx <- fnInitMatrix( rownames(lDis[[sMx]]), colnames(lDis[[sMx]]) )
        for( rr in 1:nrow(mMx)) mMx[rr,] <- as.numeric(lDis[[sMx]][rr,])
        lDis[[sMx]] <- mMx
      }
    }

    # find all risk factors by category, add each list to lDis
    lDis <- fnFindRiskFactors( lDis )
    lDis

    # convert dropped risk factors text into vector of string
    if( !is.null(lDis[["droppedRiskFactors"]])) {
      sDropped <- lDis[["droppedRiskFactors"]]
      vsDropped <- unlist(strsplit(sDropped, split="[ \t]+"))
      vbx <- !vsDropped %in% lDis$vsuAllRF
      # validate
      if( any(vbx)) stop( "Disease model 'droppedRiskFactors' specifies dropping non-existent variable(s): ", paste( vsDropped[vbx], collapse=", "))
      lDis[["droppedRiskFactors"]] <- vsDropped
    }


    # apply dropped risk factors
    if( !is.null(lDis[["droppedRiskFactors"]])) {

      # drop variables from categoricalRiskFactors
      if( !is.null(lDis[["categoricalRiskFactors"]])) {
        dfCat <- lDis[["categoricalRiskFactors"]]
        vbx <- !dfCat$riskFactor %in% lDis[["droppedRiskFactors"]]
        dfCat <- dfCat[vbx,,drop=F]
        if(!nrow(dfCat)) dfCat <- NULL
        lDis[["categoricalRiskFactors"]] <- dfCat
      }
      lDis

      # drop variables from cov matrices
      vsTotAndParts <- c("Tot",vsPartition)
      for( sDropped in lDis[["droppedRiskFactors"]] ) for( sCov in vsTotAndParts ){
        sMx <- paste0("covMatrix",sCov)
        if(is.null(lDis[[sMx]])) next
        vbxR <- !rownames(lDis[[sMx]]) %in% sDropped
        vbxC <- !colnames(lDis[[sMx]]) %in% sDropped
        lDis[[sMx]] <- lDis[[sMx]][vbxR,vbxC,drop=F]
      }
      lDis

      # update risk factors vector
      vbx <- !lDis$vsuAllRF %in% lDis[["droppedRiskFactors"]]
      lDis$vsuAllRF <- lDis$vsuAllRF[vbx]
      lDis$vsuAllRF

      # update risk factors by category, add each list to lDis
      lDis <- fnFindRiskFactors( lDis )
    }
    lDis


    if( "categoricalRiskFactors" %in% names(lDis)) {

      lifetimeRisk <- lDis$lifetimeRisk

      dfCatRf <- lDis$categoricalRiskFactors
      vsuRiskFactor <- unique(dfCatRf$riskFactor)
      bFirst <- T
      for( sRiskFactor in vsuRiskFactor) {

        vbx <- dfCatRf$riskFactor == sRiskFactor

        # validate stratum ids are unique
        if(max(table(dfCatRf$value[vbx]))!=1) stop("Problem with CategoricalRiskFactor ", sRiskFactor, " : ", "Duplicated stratum id")

        # calc per strata lifetime risks
        vStrataFreq <- dfCatRf$freq[vbx]
        vRelRisk <- dfCatRf$relRisk[vbx]
        lRF <- NA; lRF <- tryCatch(
          fnValidateCategoricalRiskFactorInputs( lifetimeRisk, vStrataFreq, vRelRisk )
          , warning=function(w) w$message, error=function(e) e$message
        ) # end of tryCatch
        if( is.character(lRF)) stop( "Problem with CategoricalRiskFactor ", sRiskFactor, " : ", lRF )

        # solve non-linear simultaneous equations
        dfRF <- lRF$dfRF
        lControl <- list(trace=0, maxit=10000)
        #fnRep("1st attempt")
        sMethod <- ifelse( nrow(dfRF)<=2,"BFGS","Nelder-Mead")
        lOptimise <- fnRunOptim( lifetimeRisk, dfRF, sMethod, lControl )
        #fnStr(lOptimise)
        if( !is.null(lOptimise$sErrMsg) & lOptimise$lOptim$sMethod!="BFGS" ) {
          # 2nd attempt
          #fnRep("2nd attempt")
          sMethod <- "BFGS"
          lOptimise <- fnRunOptim( lifetimeRisk, dfRF, sMethod, lControl )
          #fnStr(lOptimise)
        }
        if( !is.null(lOptimise$sErrMsg) ) stop("Couldn't well estimate liability distributions for strata of categorical risk factor ", sRiskFactor, " because - ", lOptimise$sErrMsg )
        dfRF <- lOptimise$dfRF
        # add missing columns to results data.frame
        if( bFirst ) {
          vsColReqd <- setdiff( colnames(dfRF), colnames(dfCatRf) )
          for( sColReqd in vsColReqd ) dfCatRf[,sColReqd] <- NA
          bFirst <- F
        }
        for( sColReqd in vsColReqd ) dfCatRf[vbx,sColReqd] <- dfRF[,sColReqd]
      }
      # calc normalised deviation from pop mean of the strata mean
      dfCatRf$liaMeanZScore <- dfCatRf$liaMean/sqrt(dfCatRf$liaVarInter)

      dfCatRf -> lDis$categoricalRiskFactors
    }
    lDis



    # create cov matrices if they don't already exist
    vsPartition <- c("Asq","Csq","Esq")
    vsTotAndParts <- c("Tot",vsPartition)
    vsCovMatrixParts <- paste0("covMatrix",vsTotAndParts)
    for( sCov in vsCovMatrixParts ) if( !sCov %in% names(lDis) ) lDis[[sCov]] <- fnInitMatrix( lDis$vsuAllRF, lDis$vsuAllRF)
    vsMCovParts <- paste0("mCov",vsTotAndParts)
    for( sCov in vsMCovParts ) lDis[[sCov]] <- fnInitMatrix( lDis$vsuAllRF, lDis$vsuAllRF)
    for( sPart in vsTotAndParts ) {
      sInput <- paste0("covMatrix",sPart)
      sMCov <- paste0("mCov",sPart)
      lDis[[sMCov]] <- fnCopyIntoMatrix(lDis[[sMCov]],lDis[[sInput]])
    }
    lDis

    for( sCov in vsCovMatrixParts ) lDis[sCov] <- NULL
    lDis

    sLia <- "liability"

    # set tot liability variance if required
    if( is.na(lDis$mCovTot[sLia,sLia]) ) lDis$mCovTot[sLia,sLia] <- 1
    if( lDis$mCovTot[sLia,sLia] != 1 ) stop("Disease liability varaiance should be set to 1")

    # copy heritability into cov matrix if required
    if( !is.na(lDis$mCovAsq[sLia,sLia]) & !is.null(lDis$heritability) ) stop("Disease heritability is specified in two places in the disease model")
    if( is.na(lDis$mCovAsq[sLia,sLia]) ) {
      if( is.null(lDis$heritability) ) stop( "Disease heritability not specified")
      lDis$mCovAsq[sLia,sLia] <- lDis$heritability * lDis$mCovTot[sLia,sLia]
    }
    # copy sharedEnv2 into cov matrix if required
    if( !is.na(lDis$mCovCsq[sLia,sLia]) & !is.null(lDis$sharedEnvironmentability) ) stop("Disease shared Environmentability is specified in two places in the disease model")
    if( is.na(lDis$mCovCsq[sLia,sLia]) ) {
      sharedEnv <- 0 # default value
      if( ! is.null(lDis$sharedEnvironmentability) ) sharedEnv <- lDis$sharedEnvironmentability
      lDis$mCovCsq[sLia,sLia] <- sharedEnv * lDis$mCovTot[sLia,sLia]
    }
    # copy uniqueEnv2 into cov matrix if required
    if( !is.na(lDis$mCovEsq[sLia,sLia]) & !is.null(lDis$uniqueEnvironmentability) ) stop("Disease unique Environmentability is specified in two places in the disease model")
    if( is.na(lDis$mCovEsq[sLia,sLia]) ) {
      uniqueEnv <- lDis$mCovTot[sLia,sLia] - (lDis$mCovAsq[sLia,sLia]+lDis$mCovCsq[sLia,sLia])  # default value
      if( ! is.null(lDis$uniqueEnvironmentability) ) uniqueEnv <- lDis$uniqueEnvironmentability
      if( uniqueEnv<0 ) stop("Disease model specifies a unique environmentability of less than zero")
      if( uniqueEnv>lDis$mCovTot[sLia,sLia] ) stop("Disease model specifies the proportion of total liability variance explained by unique environmentability exceeding 1")
      lDis$mCovEsq[sLia,sLia] <- uniqueEnv
    }
    lDis

    if( !is.null(lDis$categoricalRiskFactors)) {
      # fill in categorical risk factor variances
      for( rr in 1:nrow(lDis$categoricalRiskFactors)){
        sRF <- lDis$categoricalRiskFactors$riskFactor[rr]
        liaVarInter <- unique(lDis$categoricalRiskFactors$liaVarInter[lDis$categoricalRiskFactors$riskFactor==sRF])
        if(length(liaVarInter) != 1 ) stop("Expected exactly one value for liaVarInter")
        lDis[["mCovTot"]]["liability", sRF] <- sqrt(liaVarInter)
      }

      # fill in categorical risk factor covariances
      vsuRF <- unique(lDis$categoricalRiskFactors$riskFactor)
      sId1 <- "liability"
      for( sPart in vsPartition ){
        sMx <- paste0("mCov",sPart)
        for( sId2 in vsuRF) lDis[[sMx]][sId1,sId2] <- lDis[["mCovTot"]][sId1,sId2] * lDis[[sMx]][sId2,sId2]
      }
    }
    lDis


    # quantitative risk factors
    sId1 <- "liability"
    for( sPart in vsPartition ){
      sMx <- paste0("mCov",sPart)
      for( sId2 in lDis$vsuQuaRF) {
        if( !is.na(lDis[[sMx]][sId1,sId2]) ) stop( "A parameter that is calculated from other parameters has a value already specified for it : ", sMx, "[", sId1,",",sId2,"]" )
        lDis[[sMx]][sId1,sId2] <- lDis[["mCovTot"]][sId1,sId2] * lDis[[sMx]][sId2,sId2]
      }
    }
    lDis


    # complete cov matrices based on symmetry
    vsCovMatrixParts <- paste0("mCov",vsPartition)
    for( sCov in vsCovMatrixParts ) lDis[[sCov]] <- fnCompleteSymmetricMatrix(lDis[[sCov]]) # fill in NAs assuming matrix is symmetric

    # complete tot cov
    mCovTot <- lDis[["mCovTot"]]
    for(sR in rownames(mCovTot)) for(sC in colnames(mCovTot)) {
      totalCurr <- 0; for( sCovMatrixParts in vsCovMatrixParts) totalCurr <- totalCurr + lDis[[sCovMatrixParts]][[sR,sC]]
      if( is.na(mCovTot[sR,sC])) mCovTot[sR,sC] <- totalCurr
    }
    lDis[["mCovTot"]] <- mCovTot


    # create array
    vsTotAndParts <- c("Tot",vsPartition)
    nofRFs <- length(lDis$vsuAllRF)
    aCov <- array( NA, dim=c(nofRFs,nofRFs,length(vsTotAndParts)),dimnames=list(lDis$vsuAllRF,lDis$vsuAllRF,vsTotAndParts))
    for( sPart in vsTotAndParts ) aCov[,,sPart] <- lDis[[paste0("mCov",sPart)]][lDis$vsuAllRF,lDis$vsuAllRF,drop=F]
    lDis$aCov <- aCov

    # validate
    vsProblem <- c( vsProblem, fnValidateDis( lDis ) )

    # discard the mCov list elements
    for( sPart in vsTotAndParts ) lDis[paste0("mCov",sPart)] <- NULL

    lDis
  }

  , warning=function(w) w$message, error=function(e) e$message
  ) # end of tryCatch

  # check whether error was thrown
  if( is.character(lDisTmp)) vsProblem <- c( vsProblem, lDisTmp) else lDis <- lDisTmp

  # add problems to lDis
  if( ! is.null(vsProblem)) lDis$sErrMsg <- vsProblem

  lDis
}

#'
#' XXXX
#' @param vsLines XXXX
#' @return mCovPerson XXXX
#' @export
fnReadLinesDis <- function( vsLines )
{
  # read definition for Dis file
  dfDisDef <- diseaseRiskPredictor::dfDialectDisFile

  lDis <- list()
  vsProblem <- NULL
  lDisTmp <- tryCatch( {

    lDis <- fnReadLinesIntoDisList( vsLines )

    # do extra validation specific to dis

    # add ons
    lDis$liaThr <- stats::qnorm( p=lDis$lifetimeRisk, mean = 0, sd = 1, lower.tail = F)

    # convert the data frames read from the file into matrices
    # label the rows of the cov matrix with the id column values then drop the id column from the matrix
    vsPartition <- c("Asq","Csq","Esq")
    for( sPart in c("Tot",vsPartition) ){
      sMx <- paste0("covMatrix",sPart)
      if(!is.null(lDis[[sMx]])){

        rownames(lDis[[sMx]]) <- lDis[[sMx]]$id
        lDis[[sMx]]$id <- NULL

        # replace by matrix of the same
        vX <- apply( lDis[[sMx]], 2, as.numeric )
        mMx <- fnInitMatrix( rownames(lDis[[sMx]]), colnames(lDis[[sMx]]) )
        for( rr in 1:nrow(mMx)) mMx[rr,] <- as.numeric(lDis[[sMx]][rr,])
        lDis[[sMx]] <- mMx
      }
    }

    # find all risk factors by category, add each list to lDis
    lDis <- fnFindRiskFactors( lDis )
    lDis

    # convert dropped risk factors text into vector of string
    if( !is.null(lDis[["droppedRiskFactors"]])) {
      sDropped <- lDis[["droppedRiskFactors"]]
      vsDropped <- unlist(strsplit(sDropped, split="[ \t]+"))
      vbx <- !vsDropped %in% lDis$vsuAllRF
      # validate
      if( any(vbx)) stop( "Disease model 'droppedRiskFactors' specifies dropping non-existent variable(s): ", paste( vsDropped[vbx], collapse=", "))
      lDis[["droppedRiskFactors"]] <- vsDropped
    }


    # apply dropped risk factors
    if( !is.null(lDis[["droppedRiskFactors"]])) {

      # drop variables from categoricalRiskFactors
      if( !is.null(lDis[["categoricalRiskFactors"]])) {
        dfCat <- lDis[["categoricalRiskFactors"]]
        vbx <- !dfCat$riskFactor %in% lDis[["droppedRiskFactors"]]
        dfCat <- dfCat[vbx,,drop=F]
        if(!nrow(dfCat)) dfCat <- NULL
        lDis[["categoricalRiskFactors"]] <- dfCat
      }
      lDis

      # drop variables from cov matrices
      vsTotAndParts <- c("Tot",vsPartition)
      for( sDropped in lDis[["droppedRiskFactors"]] ) for( sCov in vsTotAndParts ){
        sMx <- paste0("covMatrix",sCov)
        if(is.null(lDis[[sMx]])) next
        vbxR <- !rownames(lDis[[sMx]]) %in% sDropped
        vbxC <- !colnames(lDis[[sMx]]) %in% sDropped
        lDis[[sMx]] <- lDis[[sMx]][vbxR,vbxC,drop=F]
      }
      lDis

      # update risk factors vector
      vbx <- !lDis$vsuAllRF %in% lDis[["droppedRiskFactors"]]
      lDis$vsuAllRF <- lDis$vsuAllRF[vbx]
      lDis$vsuAllRF

      # update risk factors by category, add each list to lDis
      lDis <- fnFindRiskFactors( lDis )
    }
    lDis


    if( "categoricalRiskFactors" %in% names(lDis)) {

      lifetimeRisk <- lDis$lifetimeRisk

      dfCatRf <- lDis$categoricalRiskFactors
      vsuRiskFactor <- unique(dfCatRf$riskFactor)
      bFirst <- T
      for( sRiskFactor in vsuRiskFactor) {
        fnStr(sRiskFactor)

        vbx <- dfCatRf$riskFactor == sRiskFactor

        # validate stratum ids are unique
        if(max(table(dfCatRf$value[vbx]))!=1) stop("Problem with CategoricalRiskFactor ", sRiskFactor, " : ", "Duplicated stratum id")

        # calc per strata lifetime risks
        vStrataFreq <- dfCatRf$freq[vbx]
        vRelRisk <- dfCatRf$relRisk[vbx]
        lRF <- NA; lRF <- tryCatch(
          fnValidateCategoricalRiskFactorInputs( lifetimeRisk, vStrataFreq, vRelRisk )
          , warning=function(w) w$message, error=function(e) e$message
        ) # end of tryCatch
        if( is.character(lRF)) stop( "Problem with CategoricalRiskFactor ", sRiskFactor, " : ", lRF )

        # solve non-linear simultaneous equations
        dfRF <- lRF$dfRF
        lControl <- list(trace=0, maxit=10000)
        fnRep("1st attempt")
        sMethod <- ifelse( nrow(dfRF)<=2,"BFGS","Nelder-Mead")
        lOptimise <- fnRunOptim( lifetimeRisk, dfRF, sMethod, lControl )
        fnStr(lOptimise)
        if( !is.null(lOptimise$sErrMsg) & lOptimise$lOptim$sMethod!="BFGS" ) {
          # 2nd attempt
          fnRep("2nd attempt")
          sMethod <- "BFGS"
          lOptimise <- fnRunOptim( lifetimeRisk, dfRF, sMethod, lControl )
          fnStr(lOptimise)
        }
        if( !is.null(lOptimise$sErrMsg) ) stop("Couldn't well estimate liability distributions for strata of categorical risk factor ", sRiskFactor, " because - ", lOptimise$sErrMsg )
        dfRF <- lOptimise$dfRF
        # add missing columns to results data.frame
        if( bFirst ) {
          vsColReqd <- setdiff( colnames(dfRF), colnames(dfCatRf) )
          for( sColReqd in vsColReqd ) dfCatRf[,sColReqd] <- NA
          bFirst <- F
        }
        for( sColReqd in vsColReqd ) dfCatRf[vbx,sColReqd] <- dfRF[,sColReqd]
      }
      # calc normalised deviation from pop mean of the strata mean
      dfCatRf$liaMeanZScore <- dfCatRf$liaMean/sqrt(dfCatRf$liaVarInter)

      dfCatRf -> lDis$categoricalRiskFactors
    }
    lDis



    # create cov matrices if they don't already exist
    vsPartition <- c("Asq","Csq","Esq")
    vsTotAndParts <- c("Tot",vsPartition)
    vsCovMatrixParts <- paste0("covMatrix",vsTotAndParts)
    for( sCov in vsCovMatrixParts ) if( !sCov %in% names(lDis) ) lDis[[sCov]] <- fnInitMatrix( lDis$vsuAllRF, lDis$vsuAllRF)
    vsMCovParts <- paste0("mCov",vsTotAndParts)
    for( sCov in vsMCovParts ) lDis[[sCov]] <- fnInitMatrix( lDis$vsuAllRF, lDis$vsuAllRF)
    for( sPart in vsTotAndParts ) {
      fnStr(sPart)
      sInput <- paste0("covMatrix",sPart)
      sMCov <- paste0("mCov",sPart)
      print(lDis[[sInput]])
      lDis[[sMCov]] <- fnCopyIntoMatrix(lDis[[sMCov]],lDis[[sInput]])
    }
    lDis

    for( sCov in vsCovMatrixParts ) lDis[sCov] <- NULL
    lDis

    sLia <- "liability"

    # set tot liability variance if required
    if( is.na(lDis$mCovTot[sLia,sLia]) ) lDis$mCovTot[sLia,sLia] <- 1
    if( lDis$mCovTot[sLia,sLia] != 1 ) stop("Disease liability varaiance should be set to 1")

    # copy heritability into cov matrix if required
    if( !is.na(lDis$mCovAsq[sLia,sLia]) & !is.null(lDis$heritability) ) stop("Disease heritability is specified in two places in the disease model")
    if( is.na(lDis$mCovAsq[sLia,sLia]) ) {
      if( is.null(lDis$heritability) ) stop( "Disease heritability not specified")
      lDis$mCovAsq[sLia,sLia] <- lDis$heritability * lDis$mCovTot[sLia,sLia]
    }
    # copy sharedEnv2 into cov matrix if required
    if( !is.na(lDis$mCovCsq[sLia,sLia]) & !is.null(lDis$sharedEnvironmentability) ) stop("Disease shared Environmentability is specified in two places in the disease model")
    if( is.na(lDis$mCovCsq[sLia,sLia]) ) {
      sharedEnv <- 0 # default value
      if( ! is.null(lDis$sharedEnvironmentability) ) sharedEnv <- lDis$sharedEnvironmentability
      lDis$mCovCsq[sLia,sLia] <- sharedEnv * lDis$mCovTot[sLia,sLia]
    }
    # copy uniqueEnv2 into cov matrix if required
    if( !is.na(lDis$mCovEsq[sLia,sLia]) & !is.null(lDis$uniqueEnvironmentability) ) stop("Disease unique Environmentability is specified in two places in the disease model")
    if( is.na(lDis$mCovEsq[sLia,sLia]) ) {
      uniqueEnv <- lDis$mCovTot[sLia,sLia] - (lDis$mCovAsq[sLia,sLia]+lDis$mCovCsq[sLia,sLia])  # default value
      if( ! is.null(lDis$uniqueEnvironmentability) ) uniqueEnv <- lDis$uniqueEnvironmentability
      if( uniqueEnv<0 ) stop("Disease model specifies a unique environmentability of less than zero")
      if( uniqueEnv>lDis$mCovTot[sLia,sLia] ) stop("Disease model specifies the proportion of total liability variance explained by unique environmentability exceeding 1")
      lDis$mCovEsq[sLia,sLia] <- uniqueEnv
    }
    lDis

    if( !is.null(lDis$categoricalRiskFactors)) {
      # fill in categorical risk factor variances
      for( rr in 1:nrow(lDis$categoricalRiskFactors)){
        sRF <- lDis$categoricalRiskFactors$riskFactor[rr]
        liaVarInter <- unique(lDis$categoricalRiskFactors$liaVarInter[lDis$categoricalRiskFactors$riskFactor==sRF])
        if(length(liaVarInter) != 1 ) stop("Expected exactly one value for liaVarInter")
        lDis[["mCovTot"]]["liability", sRF] <- sqrt(liaVarInter)
      }

      # fill in categorical risk factor covariances
      vsuRF <- unique(lDis$categoricalRiskFactors$riskFactor)
      sId1 <- "liability"
      for( sPart in vsPartition ){
        sMx <- paste0("mCov",sPart)
        for( sId2 in vsuRF) lDis[[sMx]][sId1,sId2] <- lDis[["mCovTot"]][sId1,sId2] * lDis[[sMx]][sId2,sId2]
      }
    }
    lDis


    # quantitative risk factors
    sId1 <- "liability"
    for( sPart in vsPartition ){
      sMx <- paste0("mCov",sPart)
      for( sId2 in lDis$vsuQuaRF) {
        if( !is.na(lDis[[sMx]][sId1,sId2]) ) stop( "A parameter that is calculated from other parameters has a value already specified for it : ", sMx, "[", sId1,",",sId2,"]" )
        lDis[[sMx]][sId1,sId2] <- lDis[["mCovTot"]][sId1,sId2] * lDis[[sMx]][sId2,sId2]
      }
    }
    lDis


    # complete cov matrices based on symmetry
    vsCovMatrixParts <- paste0("mCov",vsPartition)
    for( sCov in vsCovMatrixParts ) lDis[[sCov]] <- fnCompleteSymmetricMatrix(lDis[[sCov]]) # fill in NAs assuming matrix is symmetric

    # complete tot cov
    mCovTot <- lDis[["mCovTot"]]
    for(sR in rownames(mCovTot)) for(sC in colnames(mCovTot)) {
      totalCurr <- 0; for( sCovMatrixParts in vsCovMatrixParts) totalCurr <- totalCurr + lDis[[sCovMatrixParts]][[sR,sC]]
      if( is.na(mCovTot[sR,sC])) mCovTot[sR,sC] <- totalCurr
    }
    lDis[["mCovTot"]] <- mCovTot


    # create array
    vsTotAndParts <- c("Tot",vsPartition)
    nofRFs <- length(lDis$vsuAllRF)
    aCov <- array( NA, dim=c(nofRFs,nofRFs,length(vsTotAndParts)),dimnames=list(lDis$vsuAllRF,lDis$vsuAllRF,vsTotAndParts))
    for( sPart in vsTotAndParts ) aCov[,,sPart] <- lDis[[paste0("mCov",sPart)]][lDis$vsuAllRF,lDis$vsuAllRF,drop=F]
    lDis$aCov <- aCov

    # validate
    vsProblem <- c( vsProblem, fnValidateDis( lDis ) )

    # discard the mCov list elements
    for( sPart in vsTotAndParts ) lDis[paste0("mCov",sPart)] <- NULL

    lDis
  }

  , warning=function(w) w$message, error=function(e) e$message
  ) # end of tryCatch

  # check whether error was thrown
  if( is.character(lDisTmp)) vsProblem <- c( vsProblem, lDisTmp) else lDis <- lDisTmp

  # add problems to lDis
  if( ! is.null(vsProblem)) lDis$sErrMsg <- vsProblem
  fnStr(lDis)

  lDis
}

# don't call library() in a package
suppressPackageStartupMessages( library(abind) ) # for adrop()

fnValidateDis <- function( lDis )
{
  vsProblem <- NULL

  vsPartition <- c("Asq","Csq","Esq")
  vsTotAndParts <- c("Tot",vsPartition) # XXXX in fact don't need to check Tot, if partitions are +ve def then so is Tot

  # check partitions positive definite
  for( sPart in vsTotAndParts ){
    # first drop zero row/cols, otherwise tiny negative eigenvalues can occur
    mCov <- abind::adrop(lDis$aCov[,,sPart,drop=F],drop=3)
    vbxR <- apply(mCov,1,function(x){ !all(x==0)} )
    vbxC <- apply(mCov,2,function(x){ !all(x==0)} )
    mCov <- mCov[vbxR,vbxC,drop=F]

    if(!nrow(mCov) & !ncol(mCov)) next # ignore zero matrices

    lEigen <- eigen(mCov) # determine if positive definite
    if(any(lEigen$values<0)) {
      sProblem <- paste(
        "The disease model is impossible. The following cov partition is not positive definite so not a valid covariance matrix: ", sPart )
      vsProblem <- c( vsProblem, sProblem )
    }
  }
  vsProblem

  # check tot cov
  mCovTot <- abind::adrop(lDis$aCov[,,"Tot",drop=F],drop=3)
  for(sR in rownames(mCovTot))for(sC in colnames(mCovTot)){
    totalCurr <- sum( lDis$aCov[sR,sC,vsPartition] )
    diffCurr <- mCovTot[sR,sC] - totalCurr
    if( abs(diffCurr) > 1e-10 ) {
      sProblem <- paste( "The tot cov matrix is not the sum of its partitions, see row= ", sR, ", col= ", sC )
      vsProblem <- c( vsProblem, sProblem )
    }
  }

  vsProblem
}
