# don't call library() in a package
suppressPackageStartupMessages( library(abind) ) # for adrop()


#' Create population cov matrix for a single person
#'
#' XXXX
#' @param aCov XXXX
#' @return mCovPerson XXXX
#' @export
fnConvertCovArrayIntoMatrix <- function( aCov )
{
  # convert array into matrix
  dfId <- expand.grid (id1=dimnames(aCov)[[1]], id3=dimnames(aCov)[[3]] )
  dfId$full <- apply( dfId,1, paste, collapse="_" )
  mCovPerson <- fnInitMatrix( dfId$full, dfId$full, fillValue=NA)

  # the diagonal
  for( sPart in dimnames(aCov)[[3]] ){
    #XXXX mCov <- aCov[,,sPart] # R is crap, 3d array reduced to 2d matrix, unless the 2d matrix is 1x1
    mCov <- abind::adrop(aCov[,,sPart,drop=F],drop=3) # R is crap, 3d array reduced to 2d matrix, unless the 2d matrix is 1x1
    rownames(mCov) <- paste( rownames(mCov), sPart, sep="_")
    colnames(mCov) <- paste( colnames(mCov), sPart, sep="_")
    mCovPerson <- fnCopyIntoMatrix(mCovPerson,mCov)
  }
  # tot v not tot
  for( sPart in dimnames(aCov)[[3]] ){
    mCov <- abind::adrop(aCov[,,sPart,drop=F],drop=3)
    rownames(mCov) <- paste( rownames(mCov), "Tot", sep="_")
    colnames(mCov) <- paste( colnames(mCov), sPart, sep="_")
    mCovPerson <- fnCopyIntoMatrix(mCovPerson,mCov)
  }
  # not tot v not tot
  for( sPart1 in dimnames(aCov)[[3]] ) for( sPart2 in dimnames(aCov)[[3]] ){
    if( sPart1=="Tot" ) next
    if( sPart2=="Tot" ) next
    if( sPart1==sPart2 ) next
    vbx1 <- grepl( paste0("_",sPart1,"$"), rownames(mCovPerson) )
    vbx2 <- grepl( paste0("_",sPart2,"$"), colnames(mCovPerson) )
    mCovPerson[vbx1,vbx2] <- 0
  }

  mCovPerson <- fnCompleteSymmetricMatrix(mCovPerson)

  mCovPerson
}




fnAssemblePedigreeCovMatrix <- function(mCovAsqPed,mCovCsqPed,mCovEsqPed,mCovTotPed)
{
  # specify prior cov within pedigree matrix
  vsHdrCov <- c(rownames(mCovAsqPed),rownames(mCovCsqPed),rownames(mCovEsqPed),rownames(mCovTotPed))
  mCovAllPed <- matrix(NA,nrow=length(vsHdrCov),ncol=length(vsHdrCov))
  rownames(mCovAllPed) <- vsHdrCov
  colnames(mCovAllPed) <- rownames(mCovAllPed)

  # specify the block diagonal
  mCovAllPed[rownames(mCovAsqPed),rownames(mCovAsqPed)] <- mCovAsqPed
  mCovAllPed[rownames(mCovCsqPed),rownames(mCovCsqPed)] <- mCovCsqPed
  mCovAllPed[rownames(mCovEsqPed),rownames(mCovEsqPed)] <- mCovEsqPed
  mCovAllPed[rownames(mCovTotPed),rownames(mCovTotPed)] <- mCovTotPed

  # specify Tot-Other
  mCovAllPed[rownames(mCovTotPed),rownames(mCovAsqPed)] <- mCovAsqPed
  mCovAllPed[rownames(mCovTotPed),rownames(mCovCsqPed)] <- mCovCsqPed
  mCovAllPed[rownames(mCovTotPed),rownames(mCovEsqPed)] <- mCovEsqPed
  mCovAllPed[rownames(mCovTotPed),rownames(mCovTotPed)] <- mCovTotPed

  # set NonTot-Other to zero
  mbx <- row(mCovAllPed)>col(mCovAllPed) & is.na(mCovAllPed)
  mCovAllPed[mbx] <- 0

  # replace any remaining NA elements with their transpose
  mCovAllPed <- fnCompleteSymmetricMatrix( mCovAllPed )

  mCovAllPed
}

#' Calc pedigree covariance matrix for trait risk factors, liability and its partitions
#'
#' XXXX
#' @param mCovPerson XXXX
#' @param dfPed XXXX
#' @return XXXX
#' @export
fnCalcPedCov <- function( mCovPerson, dfPed )
{
  # calc correlation for each partition of (co)variance

  # calc pedigree correlation for Asq
  oPedigree <- fnConstructPedigree( dfPed )
  mCorAsq <- 2*kinship2::kinship(oPedigree)
  mCorAsq

  # calc pedigree correlation for Esq
  mCorEsq <- fnInitMatrix( dfPed$id, dfPed$id, fillValue=0)
  diag(mCorEsq) <- 1
  mCorEsq

  # calc pedigree correlation for Csq
  mCorCsq <- mCorEsq
  for(pp in 1:nrow(dfPed)) {
    if( 0==dfPed$fatherid[pp] | 0==dfPed$motherid[pp] ) next # founders can't have common (i.e.sib) env
    mCorCsq[pp,(dfPed$fatherid == dfPed$fatherid[pp] & dfPed$motherid == dfPed$motherid[pp])] <- 1
  }
  mCorCsq


  # specify population pedigree cov matrix, per partition and in total

  vbx <- NA; vbx <- grepl("_Asq",rownames(mCovPerson))
  mCovAsqPed <- kronecker(mCovPerson[vbx,vbx,drop=F],mCorAsq)
  rownames(mCovAsqPed) <- t(outer( rownames(mCovPerson[vbx,vbx,drop=F]), rownames(mCorAsq), FUN=paste, sep="_"))
  colnames(mCovAsqPed) <- rownames(mCovAsqPed)

  vbx <- NA; vbx <- grepl("_Csq",rownames(mCovPerson))
  mCovCsqPed <- kronecker(mCovPerson[vbx,vbx,drop=F],mCorCsq)
  rownames(mCovCsqPed) <- t(outer( rownames(mCovPerson[vbx,vbx,drop=F]), rownames(mCorCsq), FUN=paste, sep="_"))
  colnames(mCovCsqPed) <- rownames(mCovCsqPed)

  vbx <- NA; vbx <- grepl("_Esq",rownames(mCovPerson))
  mCovEsqPed <- kronecker(mCovPerson[vbx,vbx,drop=F],mCorEsq)
  rownames(mCovEsqPed) <- t(outer( rownames(mCovPerson[vbx,vbx,drop=F]), rownames(mCorEsq), FUN=paste, sep="_"))
  colnames(mCovEsqPed) <- rownames(mCovEsqPed)

  # calc total liability pedigree cov
  mCovTotPed <- mCovAsqPed + mCovCsqPed + mCovEsqPed
  rownames(mCovTotPed) <- sub("_Asq_","_Tot_",rownames(mCovTotPed))
  colnames(mCovTotPed) <- rownames(mCovTotPed)


  # assemble per partition and total pedigree covs
  fnAssemblePedigreeCovMatrix( mCovAsqPed, mCovCsqPed, mCovEsqPed, mCovTotPed)
}



#PPPP
# calculates posterior distribution via Pearson Aitken selection formula
# C = all variables comprising 2 parts A and B
# A = variables which have undergone selection
# B = variables of C that are not in A
fnPearsonAitken <- function( priorMeanC,priorCovC, postMeanA,postCovA )
{
  if( !is.matrix(priorMeanC) ) stop( "expected matrix for prior mean")
  if( !is.matrix(postMeanA) ) stop( "expected matrix for posterior mean")

  vsC <- rownames(priorMeanC)
  vsA <- rownames(postMeanA)
  if( is.null(vsC) ) stop( "expected rownames for prior mean")
  if( is.null(vsA) ) stop( "expected rownames for posterior mean")

  if( !all(vsA %in% vsC)) stop( "Expected A to be a subset of C \n")

  # calc B set
  vsB <- vsC[ !vsC %in% vsA]

  # validate
  if( nrow(priorCovC) != ncol(priorCovC) ) stop( "Expected square matrix cov C \n")
  if( !all(rownames(priorMeanC) == rownames(priorCovC) )) stop( "Expected rownames to match for C \n")
  if( nrow(postCovA) != ncol(postCovA) ) stop( "Expected square matrix cov A \n")
  if( !all(rownames(postMeanA) == rownames(postCovA) )) stop( "Expected rownames to match for A \n")

  # calc post cov B
  priorInvCovA <- solve(priorCovC[vsA,vsA,drop=F])
  priorCovBA <- priorCovC[vsB,vsA,drop=F]
  mRegCoeff <- priorCovBA %*% priorInvCovA
  colnames(mRegCoeff) <- vsA

  postCovB <- priorCovC[vsB,vsB,drop=F] - priorCovC[vsB,vsA,drop=F] %*%
    ( priorInvCovA - priorInvCovA %*% postCovA %*% priorInvCovA ) %*% priorCovC[vsA,vsB,drop=F]
  # calc post cov AB
  postCovAB <- postCovA %*% priorInvCovA %*% priorCovC[vsA,vsB,drop=F]
  # calc post cov C
  postCovC <- rbind( cbind( postCovA, postCovAB), cbind( t(postCovAB), postCovB))
  postCovC <- postCovC[vsC,vsC,drop=F]

  # calc post mean B
  postMeanB <- priorMeanC[vsB,,drop=F] + priorCovC[vsB,vsA,drop=F] %*% priorInvCovA %*% ( postMeanA - priorMeanC[vsA,,drop=F] )
  # calc post mean C
  postMeanC <- rbind(postMeanA, postMeanB, drop=F )
  postMeanC <- postMeanC[vsC,,drop=F]

  list( postMeanC=postMeanC, postCovC=postCovC, vsB=vsB, mRegCoeff=mRegCoeff )
}



#' XXXX
#'
#' XXXX
#' @param vMeanAllPed XXXX
#' @param mCovAllPed XXXX
#' @param vObsRFValues XXXX
#' @return XXXX
#' @export
fnCalcPedDistributionPosteriorToObsRiskFactors <- function( vMeanAllPed, mCovAllPed, vObsRFValues )
{
  mPriorCov <- mCovAllPed
  mPriorMean <- as.matrix(vMeanAllPed,ncol=1)

  if(!length(vObsRFValues)) return(list( mPostMean=mPriorMean, mPostCov=mPriorCov ))

  #### condition on observed risk factors
  vsA <- names(vObsRFValues)
  postMeanA <- vObsRFValues
  postMeanA
  postCovA <- mPriorCov[vsA,vsA,drop=F]; postCovA[,] <- 0
  postCovA
  lPA <- fnPearsonAitken( mPriorMean, mPriorCov, as.matrix(postMeanA,ncol=1), postCovA )
  # get posterior distribution dropping all variables with no uncertainty
  vbxHasVariance <- diag(lPA$postCovC)!=0
  mPostCov <- lPA$postCovC[vbxHasVariance,vbxHasVariance]
  mPostMean <- (lPA$postMeanC)[vbxHasVariance,1,drop=F]

  list( mPostMean=mPostMean, mPostCov=mPostCov )
}



#' Collate observations on risk factors
#'
#' XXXX
#' @param lDis XXXX
#' @param dfPed XXXX
#' @return XXXX
#' @export
fnCollateObservedRiskFactors <- function( lDis, dfPed )
{
  # categorical risk factors
  vsuCatRf <- unique(lDis$categoricalRiskFactors$riskFactor)
  vsuCatRf <- vsuCatRf[ !vsuCatRf %in% lDis$droppedRiskFactors ]

  # quantitative risk factors
  vsRF <- rownames(lDis$aCov)
  vbx <- !vsRF %in% c(vsuCatRf,"liability")
  vsuQuantRf <- vsRF[vbx]

  # validate
  if( length(intersect(vsuCatRf, vsuQuantRf))) stop( "The following risk factors are defined as both categorical and quantitative: ", paste( intersect(vsuCatRf, vsuQuantRf), collapse=", "))

  # validate risk factors are in pedigree
  vsMissingFromPed <- setdiff( c(vsuCatRf, vsuQuantRf), colnames(dfPed) )
  if( length(vsMissingFromPed) ) stop("The pedigree file is missing variables defined as risk factors in the disease model. These are - ", paste( vsMissingFromPed, collapse=" " ) )

  dfRFObs <- dfPed[ , vsuCatRf, drop=F ]
  dfRFObs
  # replace any categorical risk factors by appropriate value to cause specified rel risk
  for(rr in 1:nrow(dfRFObs)) for(sCatRf in vsuCatRf) {
    sValue <- dfRFObs[rr,sCatRf]
    if( is.na(sValue) ) next

    vbx <- lDis$categoricalRiskFactor$riskFactor == sCatRf & lDis$categoricalRiskFactor$value == sValue
    nofMatches <- sum(vbx)
    # validate
    if(!nofMatches) stop("In the pedigree file a categorical risk factor has been given a value that is not defined in the disease model", "Categorical risk factor ='", sCatRf, "'")
    if(nofMatches>1) stop("In the pedigree file a categorical risk factor has been given a value that is multiply defined in the disease model", "Categorical risk factor ='", sCatRf, "'")

    # replace categorical value with value to specified rel risk
    dfRFObs[rr,sCatRf] <- lDis$categoricalRiskFactors$liaMeanZScore[vbx]
  }
  rownames(dfRFObs) <- paste( "Tot", dfPed$id, sep="_")
  dfRFObs

  # convert into named vector
  vObsRFValues <- as.vector(as.matrix(dfRFObs))
  names(vObsRFValues) <- as.vector(t(outer( colnames(dfRFObs), rownames(dfRFObs), paste, sep="_")))
  vObsRFValues


  # quantitative risk factors
  dfRFObs <- dfPed[ , vsuQuantRf, drop=F ]
  dfRFObs

  # check pedigree info for quantitative risk factors is numeric
  vbIsNotNumeric <- ! apply( dfRFObs,2,is.numeric )
  vbxAllNA <- apply(dfRFObs,2, function(x){ all(is.na(x)) } )
  vsNotNumeric <- colnames(dfRFObs)[vbIsNotNumeric & !vbxAllNA]
  if( length(vsNotNumeric) ) stop( "Risk factors specified as quantitative in the disease model are represented by non-numeric values in the pedigree information. The problem risk factors are: ", paste( vsNotNumeric, collapse=", "))

  rownames(dfRFObs) <- paste( "Tot", dfPed$id, sep="_")
  # convert into named vector
  vObsRFValuesQuant <- as.vector(as.matrix(dfRFObs))
  names(vObsRFValuesQuant) <- as.vector(t(outer( colnames(dfRFObs), rownames(dfRFObs), paste, sep="_")))
  vObsRFValuesQuant

  vObsRFValues <- c( vObsRFValues, vObsRFValuesQuant)

  # discard unobserved variables
  vObsRFValues <- vObsRFValues[!is.na(vObsRFValues)]
  vObsRFValues
}



# calculates parameters needed for calculating each variable's distribution conditional on a point value for the remaining variables
# for each variable, the parameters calculated are the reg coeff vector and posterior std dev
fnCalcConditioningOnRestOfVariables <- function( mCov )
{
  vsCovVars <- rownames(mCov)
  mMean <- matrix(0,length(vsCovVars),1); rownames(mMean) <- vsCovVars
  mCov0 <- mCov; mCov0[] <- 0; mCov0
  lmRegCoeff <- list()
  vPostStdDevB <- rep( NA, nrow(mCov)); names(vPostStdDevB) <- vsCovVars
  for( bb in 1:length(vsCovVars) ) {
    sB <- vsCovVars[bb]
    vsA <- vsCovVars[-bb,drop=F]
    mMeanA <- mMean[vsA,1,drop=F]
    mCovA <- mCov0[vsA,vsA,drop=F]
    lPA <- fnPearsonAitken( mMean, mCov, mMeanA, mCovA )

    vPostStdDevB[sB] <- sqrt(lPA$postCovC[lPA$vsB,lPA$vsB])
    lmRegCoeff[[sB]] <- lPA$mRegCoeff
  }
  list(lmRegCoeff=lmRegCoeff, vPostStdDevB=vPostStdDevB)
}

#' XXXX
#'
#' XXXX
#' @param dfPed XXXX
#' @param mPostCov XXXX
#' @param mPostMean XXXX
#' @param nofBurnIn XXXX
#' @param nofIter XXXX
#' @return XXXX
#' @export
fnSampleConditionedLiaDistribution <- function( dfPed, mPostCov, mPostMean, nofBurnIn=100, nofIter=9000 )
{
  lConditioningOnRest <- fnCalcConditioningOnRestOfVariables( mPostCov )
  lmRegCoeff <- lConditioningOnRest$lmRegCoeff
  vPostStdDevB <- lConditioningOnRest$vPostStdDevB

  bUseCppVersion <- T
  #set.seed(42)
  fnGibbsSampler( dfPed, mPostMean, lmRegCoeff, vPostStdDevB, nofBurnIn, nofIter, bUseCppVersion )
}



# GGGG
# don't call library() in a package
suppressPackageStartupMessages( library(Rcpp))
suppressPackageStartupMessages( library(truncnorm))

fnGibbsSampler <- function( dfPed, mPriorMean, lmRegCoeff, vPostStdDev, nofBurnIn, nofIter, useCpp = TRUE, iMaxAttempts=10000000 )
{
  #fnRep("fnGibbsSampler()")
  #fnStr(useCpp)
  # validate input types
  if( !is.matrix(mPriorMean)) stop("Expected a matrix for - mPriorMean")
  if( !is.list(lmRegCoeff)) stop("Expected a list for - lmRegCoeff")
  if( !is.vector(vPostStdDev)) stop("Expected a vector for - vPostStdDev")

  # validate dfPed
  vsRequiredPedCol <- c("id","affected","expressedProportionOfLifetimeRisk","thr")
  vbx <- ! vsRequiredPedCol %in% colnames(dfPed)
  if( any(vbx)) stop( "dfPed missing the following columns: ", paste( vsRequiredPedCol[vbx], collapse=", ") )

  nofVariables <- length(as.vector(mPriorMean))
  vsVariables <- rownames(mPriorMean)

  # validate compatible inputs
  if(nofVariables!=length(lmRegCoeff)|nofVariables!=length(vPostStdDev)) stop("Inconsistent input parameters for function fnGibbsSampler() - differing lengths")
  if(any(vsVariables!=names(lmRegCoeff))|any(vsVariables!=names(vPostStdDev))) stop("Inconsistent input parameters for function fnGibbsSampler() - differing names")

  # validate
  if(1<sum(0==vPostStdDev)) stop("Cannot use Gibbs sampler when conditional std dev is zero for any pair or more of variables")

  # parse variable names
  lvsVariables <- strsplit(vsVariables,"_")
  msParseVariables <- t(matrix(unlist(lvsVariables),ncol=length(lvsVariables)))

  # index liability variables and non liability variables
  vbxLia <- "liability" == msParseVariables[,1]
  # get person ids
  vPerson <- msParseVariables[,3]

  #fnStr(dfPed); fnStr(mPriorMean); fnStr(lmRegCoeff); fnStr(vPostStdDev); fnStr(nofBurnIn); fnStr(nofIter)
  #fnStr(nofVariables); fnStr(vsVariables); fnStr(vPerson); fnStr(vbxLia)

  # specify starting point in parameter space for Gibbs Sampler
  if(F){
    vDraw <- mPriorMean[,1]
    meanLiaUna <- truncnorm::etruncnorm(a=-Inf, b=unique(dfPed$thr), mean=0, sd=1 )
    meanLiaAff <- truncnorm::etruncnorm(a=unique(dfPed$thr), b=Inf, mean=0, sd=1 )
    vDraw <- mPriorMean[,1]
    vbx <- dfPed$affected==0; vbx[is.na(vbx)] <- F
    vDraw[vbx] <- meanLiaUna
    vbx <- dfPed$affected==1; vbx[is.na(vbx)] <- F
    vDraw[T] <- meanLiaAff
  }
  vDraw <- mPriorMean[,1]
  vDraw[] <- dfPed$thr
  vDraw

  if (useCpp) {
    #Rcpp::sourceCpp("src/rcppFunctions.cpp")
    mDraws <- cppFnGibbsSampler(dfPed, mPriorMean, lmRegCoeff, vPostStdDev, vDraw, nofBurnIn, nofIter, nofVariables, vsVariables, vPerson, vbxLia, iMaxAttempts)

  } else {

    # init draws repository
    mDraws <- matrix( NA, nrow=nofIter, ncol=nofVariables ); colnames(mDraws) <- vsVariables

    for( ii in -nofBurnIn:nofIter) {
      if( ii %% 1000 == 0 ) { cat("INFO: iteration", ii, "of", nofIter, "\n" ); utils::flush.console() }

      for( vv in 1:length(vsVariables) ) {

        sB <- vsVariables[vv]
        vsA <- vsVariables[-vv,drop=F] # the set of variables to be conditioned on

        # calc posterior distribution
        meanB <- mPriorMean[sB,1] + lmRegCoeff[[sB]] %*% (as.matrix(vDraw[vsA],ncol=1)-mPriorMean[vsA,1,drop=F])
        stdDevB <- vPostStdDev[sB]

        if( !vbxLia[vv] ){
          # draw sample
          vDraw[sB] <- stats::rnorm( n=1, mean=meanB, sd=stdDevB)

        } else {
          # draw sample of total liability conditioned on affection status and expressed proportion of lifetime risk
          # find index for current person in the ped
          sPerson <- vPerson[vv]
          pp <- which(sPerson == dfPed$id)
          if(1!=length(pp)) stop("Expected exactly one person in pedigree with id = '", sPerson, "'")

          # draw a sample that matches the person's observed affection status
          for(ss in 1:iMaxAttempts) {

            # draw sample
            liaB <- stats::rnorm( n=1, mean=meanB, sd=stdDevB)

            # draw affection status based on total liability
            bAffB <- liaB > dfPed$thr[pp]
            if(bAffB) bAffB <- stats::runif(1) < dfPed$expressedProportionOfLifetimeRisk[pp]

            # check whether drawn affection status matches observed
            if( is.na(dfPed$affected[pp]) ) break
            if( dfPed$affected[pp] == 0 && !bAffB ) break
            if( dfPed$affected[pp] == 1 &&  bAffB ) break
          }
          # report problem if maximum number of attempts exceeded
          if( ss == iMaxAttempts ) stop( "The affection status pattern in the pedigree seems very unlikely given the disease model. Unable to draw a liability sample for person ", sPerson, " despite ", iMaxAttempts, " attempts" )

          # record draw from person's liability distribution
          vDraw[sB] <- liaB
        }
      }

      # record draw from pedigree's liability distribution
      if(ii>0) mDraws[ii,] <- vDraw
    }

  }

  return(mDraws)
}



#' Calculates risks based on a sample liability dataset
#'
#' XXXX
#' @param dfPed XXXX
#' @param mTotLiaSample XXXX
#' @return XXXX
#' @export
fnEstimateRiskFromTotLiaSample <- function( dfPed, mTotLiaSample )
{
  dfPed$risk <- NA
  for(rr in 1:nrow(dfPed)) {
    sPerson <- dfPed$id[rr]
    sColTotLia <- paste("liability","Tot", sPerson, sep="_")
    vbx <- colnames(mTotLiaSample) == sColTotLia
    if( !sum(vbx) ) next
    if( sum(vbx)>1 ) stop("Expected exactly 1 column in mTotLiaSample called ", sColTotLia )
    dfPed$risk[rr] <- mean( mTotLiaSample[,sColTotLia] > dfPed$thr[rr] )
  }
  dfPed
}



#' Calculates expressed proportion of lifetime risk for pedigree members
#'
#' If disease model congenital then EPLR=1 for all. Otherwise
#' Attaches an age of onset curve to each person according to their personal attributes.
#' Interpolates each persons EPLR from their age and age of onset curve
#' @param lDis disease model
#' @param dfPed pedigree
#' @return pedigree with expressedProportionOfLifetimeRisk
#' @export
fnCalcExprPropLifetimeRisk <- function(lDis, dfPed)
{
  #fnRep("fnCalcExprPropLifetimeRisk")

  # calculate pedigree members' expressed proportion of lifetime risk
  sParam <- "ageOfOnset"
  if( ! sParam %in% names(lDis) ) {
    # disease is congenital
    dfPed$expressedProportionOfLifetimeRisk <- 1
    return(dfPed)
  }

  # validate pedigree age
  if( ! "age" %in% colnames(dfPed) ) stop( "Can't apply age of onset information. No 'age' variable is specified in the pedigree data" )
  if( ! is.numeric(dfPed$age) ) stop( "Can't apply age of onset information. The 'age' variable specified in the pedigree data is not numeric" )
  if( any(is.na(dfPed$age) | is.nan(dfPed$age) | is.infinite(dfPed$age)) ) stop( "Can't apply age of onset information. The 'age' variable specified in the pedigree data is not a good number for at least one pedigree member" )

  dfPed$age
  dfAOO <- lDis[[sParam]]
  # get AOO stratum defining variables
  vbxAOOStratumDefVariables <- ! colnames(dfAOO) %in% c( "age","expressedProportionOfLifetimeRisk")
  vsStratumDefVariables <- colnames(dfAOO)[vbxAOOStratumDefVariables]
  # validate A00 statum defining variables are in ped
  if(length(vsStratumDefVariables)){
    vbx <- !vsStratumDefVariables %in% colnames(dfPed)
    if(any(vbx)) stop( "The pedigree does not contain the following age of onset stratum defining variables: ", paste( vsStratumDefVariables[vbx], collapse=", "))
  }

  dfPed$expressedProportionOfLifetimeRisk <- NA
  for( pp in 1:nrow(dfPed)) {
    #fnStr(pp)
    #fnStr(dfPed[pp,])
    # find AOO for current person
    vbx <- NA
    if(!length(vsStratumDefVariables)){
      vbx <- T # a single AOO stratum
    } else { # obtain approp AOO stratum
      vbx <- apply( dfAOO[,vsStratumDefVariables,drop=F], 1, function(x){ all(x == dfPed[pp,vsStratumDefVariables] )})
      # validate
      if( !sum(vbx) ) stop( "An age of onset curve is not defined for pedigree member: ", "familyId=", dfPed$famid, ", ", "personId=", dfPed$id )
    }
    dfAOOcurr <- dfAOO[vbx,,drop=F]
    #fnStr(dfAOOcurr)
    # validate
    if(length(unique(dfAOOcurr$age))<2) stop( "Age of onset curve is inadequate (not enough ages) for pedigree member: ", "familyId=", dfPed$famid, ", ", "personId=", dfPed$id )

    # find closest pair of AOO ages spanning person's age
    dfAOOcurr <- dfAOOcurr[ order(dfAOOcurr$age),]
    vDiffAge <- dfPed[pp,"age"] - dfAOOcurr$age
    vbDiffAgePositive <- vDiffAge >= 0
    vPos <- which(vbDiffAgePositive)
    vNeg <- which(!vbDiffAgePositive)
    iStart <- NA
    if(length(vPos)) iStart <- max(vPos)
    if(is.na(iStart)) iStart <- min(vNeg)
    iStop <- iStart + 1
    vixNearest <- iStart:iStop
    if(iStop>nrow(dfAOOcurr)) vixNearest <- vixNearest - 1
    #fnStr(dfAOOcurr[vixNearest,])

    # interpolate expressedProportionOfLifetimeRisk for person
    reg <- stats::lm(formula=expressedProportionOfLifetimeRisk ~ age, data=dfAOOcurr[vixNearest,])
    expressedProportionOfLifetimeRisk <- stats::coef(reg)["(Intercept)"] + stats::coef(reg)["age"] * dfPed[pp,"age"]
    #fnStr(expressedProportionOfLifetimeRisk)
    # correct for bounds
    if(expressedProportionOfLifetimeRisk>1) expressedProportionOfLifetimeRisk <- 1
    if(expressedProportionOfLifetimeRisk<0) expressedProportionOfLifetimeRisk <- 0
    # assign
    dfPed$expressedProportionOfLifetimeRisk[pp] <- expressedProportionOfLifetimeRisk
    #fnStr(dfPed)
  }
  dfPed
}



# don't call library() in a package
suppressPackageStartupMessages(library(graphics))
suppressPackageStartupMessages(library(plotrix))

#' XXXX
#'
#' XXXX
#' @param dfPedRisk pedigree with disease risks
#' @return NULL
#' @export
fnBarPlotRisk <- function( dfPedRisk )
{
  vRisk <- dfPedRisk$risk
  vId <- dfPedRisk$id
  if(!is.vector(vRisk)|!is.vector(vId)) stop('Expected vectors')
  if(length(vRisk)!=length(vId)) stop('Expected vectors of same length')
  if(!length(vId)) stop('Zero length')

  secondHighestRiskMaxForGap <- 0.4
  vbxRisk1 <- !( vRisk!=1 | is.na(vRisk) )
  secondHighestRisk <- max(vRisk[!vbxRisk1],na.rm=T)
  bUseGapPlot <- secondHighestRisk<secondHighestRiskMaxForGap & sum(vbxRisk1)>0

  sMain <- "Disease Risk in Pedigree"
  sX <- "Person Id"
  sY <- "Risk"
  if(!bUseGapPlot){
    graphics::barplot( vRisk, names.arg=vId,
             main=sMain, xlab=sX, ylab=sY )
  } else {
    vGap <- c( secondHighestRiskMaxForGap+0.05, 0.7 )
    plotrix::gap.barplot( y=vRisk, xaxlab=vId,
                 main=sMain, xlab=sX, ylab=sY,
                 gap=vGap, col=rep("grey",nrow(dfPedRisk)) )
  }
}




#' XXXX
#'
#' XXXX
#' @param dfPed pedigree
#' @param lDis disease model
#' @param bVerbose boolean for reporting plots
#' @return list containing liability draws, risk estimates, etc
#' @export
fnPrepareForRiskPrediction <- function( dfPed, lDis, bVerbose=F )
{
  lOut <- tryCatch({

    lOut <- list( dfPed=dfPed, lDis=lDis )

    # pedigree
    nofPersons <- nrow(dfPed)
    vsPerson <- dfPed$id

    # create population cov matrix for a single person of trait liability and risk factors for total liability and its partitions
    #fnRep("create population cov matrix for a single person")
    mCovPerson <- fnConvertCovArrayIntoMatrix( lDis$aCov )
    mCovPerson
    if(bVerbose) fnHeatMapCovMatrix(mCovPerson)

    #### calc pedigree distribution for trait and risk factors, liability and its partitions
    #fnRep("calc distribution for pedigree")

    # calc pedigree covariance matrix for trait risk factors, liability and its partitions
    mCovAllPed <- fnCalcPedCov( mCovPerson, dfPed )
    # drop variables with zero variance
    vbxR <- apply(mCovAllPed,1,function(x){ !all(x==0)} )
    vbxC <- apply(mCovAllPed,2,function(x){ !all(x==0)} )
    mCovAllPed <- mCovAllPed[vbxR,vbxC,drop=F]

    # specify population pedigree mean vector
    vMeanAllPed <- rep(0,nrow(mCovAllPed)); names(vMeanAllPed) <- rownames(mCovAllPed)

    if(bVerbose) fnHeatMapCovMatrix(mCovAllPed)
    #fnStr(mCovAllPed)
    #fnStr(vMeanAllPed)
    lOut$mPopCov <- mCovAllPed
    lOut$mPopMean <- as.matrix(vMeanAllPed)

    #### collate observations on risk factors
    #fnRep("collate observations on risk factors")
    vObsRFValues <- fnCollateObservedRiskFactors( lDis, dfPed )
    #fnStr(vObsRFValues)
    #print(vObsRFValues)



    # condition pedigree distribution on the observed risk factors
    #fnRep("condition pedigree distribution on the observed risk factors")
    lPedDistributionPosterior <- fnCalcPedDistributionPosteriorToObsRiskFactors( vMeanAllPed, mCovAllPed, vObsRFValues )
    mPostMean <- lPedDistributionPosterior$mPostMean
    mPostCov <- lPedDistributionPosterior$mPostCov
    if(bVerbose) fnHeatMapCovMatrix(mPostCov)
    #fnStr(mPostCov)
    #fnStr(mPostMean)

    # drop variables other than total variables
    # this is to ensure a non-singular cov matrix, if tot and its partitions are included the matrix will be singular
    vbxTot <- grepl("_Tot",rownames(mPostMean))
    mPostMean <- mPostMean[vbxTot,1,drop=F]
    mPostCov <- mPostCov[vbxTot,vbxTot]
    if(bVerbose) fnHeatMapCovMatrix(mPostCov)

    # drop variables other than total liability variables
    vbxTot <- grepl("^liability_Tot",rownames(mPostMean))
    mPostMean <- mPostMean[vbxTot,1,drop=F]
    mPostCov <- mPostCov[vbxTot,vbxTot]
    #fnRep("The mean and cov of the pedigree distribution conditioned on the observed risk factors follow")
    #fnRep("Mean =")
    mPostMean
    #fnRep("Cov =")
    mPostCov
    lOut$mPriorCov <- mPostCov
    lOut$mPriorMean <- mPostMean



    # copy per person disease model parameters from the disease model to the pedigree info
    #fnRep("copy per person disease model parameters from the disease model to the pedigree info")
    sParam <- "lifetimeRisk"
    if( ! sParam %in% names(lDis) ) stop("Disease model missing expected parameter - ", sParam)
    dfPed[,sParam] <- lDis$lifetimeRisk
    #dfPed$thr <- stats::qnorm(1-dfPed$lifetimeRisk)
    dfPed$thr <- lDis$liaThr #XXXX replace $thr with $liaThr

    dfPed <- fnCalcExprPropLifetimeRisk(lDis, dfPed)
    #fnStr(dfPed)

    lOut$dfPed <- dfPed

    lOut
  }

  , warning=function(w) {
    lOut$sWrnMsg <- w$message
    lOut$w <- w; fnStr(w)
    lOut
  }
  , error=function(e) {
    lOut$sErrMsg <- e$message
    e$traceback <- traceback()
    lOut$e <- e; fnStr(e)
    lOut
  }
  ) # end of tryCatch

  lOut
}




#' Predicts risk in disease pedigree
#'
#' A Gibbs sampler is used to draw samples from the pedigree's posterior liability distribution.
#' From this set of draws risk and n year risk is predicted for each pedigree member.
#' @param lPedDis disease model and pedigree info
#' @param nofBurnIn number of initial iterations of Gibbs sampler to be ignored
#' @param nofDraws number of draws from posterior pedigree liability distribution
#' @param nofYears number of years for n year risk calculation
#' @param bVerbose boolean for reporting plots
#' @return list containing liability draws, risk estimates, etc
#' @export
fnPredictRisk <- function( lPedDis, nofBurnIn=100, nofDraws=30000, nofYears=5, bVerbose=F )
{
  lRes <- tryCatch({

    lRes <- lPedDis

    dfPed <- lPedDis$dfPed
    lDis <- lPedDis$lDis
    mPriorCov <- lPedDis$mPriorCov
    mPriorMean <- lPedDis$mPriorMean

    #### obtain samples (via Gibbs Sampling) of the pedigree distribution conditioned on affection status and expressed proportion of lifetime risk
    mTotLiaSample <- NA; mTotLiaSample <-
      fnSampleConditionedLiaDistribution( dfPed, mPriorCov, mPriorMean, nofBurnIn=nofBurnIn, nofIter=nofDraws )
    lRes$mTotLiaSample <- mTotLiaSample

    # posterior pedigree distribution, compare to inputs
    lRes$mPostCov <- stats::cov(mTotLiaSample)
    lRes$mPostMean <- as.matrix(apply(mTotLiaSample,2,mean))

    # calc autocorrelation
    iLagMax <- 20
    nofIter <- nrow(mTotLiaSample)
    mAutoCor <- matrix( NA, nrow=iLagMax, ncol=ncol(mTotLiaSample) ); colnames(mAutoCor) <- colnames(mTotLiaSample)
    for(iLag in 1:iLagMax ) mAutoCor[iLag,] <- diag(stats::cor(mTotLiaSample[-(1:iLag),], mTotLiaSample[-((1+nofIter-iLag):nofIter),]))
    lRes$mAutoCor <- mAutoCor

    # risk predictions for the pedigree members
    dfPedRisk <- fnEstimateRiskFromTotLiaSample( dfPed, mTotLiaSample )

    # n year risk predictions for the pedigree members
    dfPedRisk <- fnCalcNYearRisk( dfPedRisk, lDis, mTotLiaSample, nofYears )

    lRes$dfPedRisk <- dfPedRisk

    lRes
  }

  , warning=function(w) {
    lRes$sWrnMsg <- w$message
    lRes$w <- w; fnStr(w)
    lRes
  }
  , error=function(e) {
    lRes$sErrMsg <- e$message
    e$traceback <- traceback()
    lRes$e <- e; fnStr(e)
    lRes
  }
  ) # end of tryCatch

  lRes
}



fnCalcNYearRisk <- function( dfPed, lDis, mTotLiaSample, nofYears=5 )
{
  dfPed$nofYears <- nofYears

  # calc expressed proportion of lifetime risk for each pedigree member at a later age
  dfPedLater <- dfPed
  dfPedLater$age <- dfPedLater$age + nofYears
  dfPedLater <- fnCalcExprPropLifetimeRisk( lDis, dfPedLater )

  # calc paramters
  dfDeltaRisk <- data.frame(
    famid = dfPedLater$famid, id = dfPedLater$id,
    f1 = 1 - dfPed$expressedProportionOfLifetimeRisk,
    f2 = 1 - dfPedLater$expressedProportionOfLifetimeRisk,
    a=NA, b=NA )
  for( rr in 1:nrow(dfDeltaRisk)){
    #XXXX NB indexed by person id so must be unique in ped
    vDraw <- mTotLiaSample[,paste0("liability_Tot_", dfPed$id[rr])]
    dfDeltaRisk$a[rr] <- sum(vDraw <  dfPed$thr[rr])
    dfDeltaRisk$b[rr] <- sum(vDraw >= dfPed$thr[rr])
  }
  dfDeltaRisk$nYearRisk <- with( dfDeltaRisk, expr= b*(1-f2/f1)/(a+b) )

  dfPed <- merge(dfPed, dfDeltaRisk, by = c("famid","id") )

  # n year risk for affected
  vbx <-  !(!dfPed$affected | is.na(dfPed$affected))
  dfPed$nYearRisk[vbx] <- 0

  # n year risk for unknown affection status
  vbx <- is.na(dfPed$affected)
  dfPed$nYearRisk[vbx] <- NA

  dfPed
}
