#' Partial valdidation of a pedigree data.frame
#'
#' Partial valdidation of a pedigree data.frame
#' @param dfPed pedigree data.frame
#' @return dfPed
#' @export
fnValidatePedigree <- function( dfPed )
{
  #fnRep("fnValidatePedigree")
  #fnStr(class( data.frame(letters)[,1]))

  # validate no duplicate column headings
  tblHeadingsPed <- table(colnames(dfPed))
  vbxHeadingsPedDuplicated <- tblHeadingsPed !=1
  if( sum(vbxHeadingsPedDuplicated)){
    sErrMsg <- paste( "Ped file contains duplicate columns headings  -", paste( names(tblHeadingsPed)[vbxHeadingsPedDuplicated], collapse=", ") )
    stop(sErrMsg)
  }

  dfValid <- rbind(
    data.frame(  col="famid",    sFnValidate=NA,                        min=NA, max=NA, bAllowNAs=F, bRequired=T )
    ,data.frame( col="id",       sFnValidate=NA,                        min=NA, max=NA, bAllowNAs=F, bRequired=T )
    ,data.frame( col="fatherid", sFnValidate=NA,                        min=NA, max=NA, bAllowNAs=F, bRequired=T )
    ,data.frame( col="motherid", sFnValidate=NA,                        min=NA, max=NA, bAllowNAs=F, bRequired=T )
    ,data.frame( col="sex",      sFnValidate="fnValidateIntegerVector", min=1, max=4,   bAllowNAs=F, bRequired=T )
    # XXXX kinship2 bug - pedigree checking allows affected=NA but plot.pedigree crashes
    #,data.frame( col="affected", sFnValidate="fnValidateIntegerVector", min=0, max=1,   bAllowNAs=F, bRequired=T )
    ,data.frame( col="affected", sFnValidate="fnValidateIntegerVector", min=0, max=1,   bAllowNAs=T, bRequired=T )
    ,data.frame( col="age",      sFnValidate="fnValidateNumericVector", min=0, max=150, bAllowNAs=T, bRequired=F )
    ,data.frame( col="deceased", sFnValidate="fnValidateIntegerVector", min=0, max=1,   bAllowNAs=F, bRequired=F )
  )
  dfValid

  # validate pedigree info contains all required columns
  vsHeadingReqdMissing <- setdiff( dfValid$col[dfValid$bRequired], colnames(dfPed))
  if( length(vsHeadingReqdMissing) ) stop( paste( "Ped file missing the following required columns -", paste( vsHeadingReqdMissing, collapse=", ") ) )

  # validate each column for which there is validation info
  #fnStr(dfPed)
  for( rr in 1:nrow(dfValid)){
    sCol <- dfValid$col[rr]
    if( !sCol %in% colnames(dfPed) ) next

    # validate for NAs
    if( (! dfValid$bAllowNAs[rr]) & any(is.na(dfPed[,sCol])) ) stop( "Problem with pedigree column '",sCol,"': ", "NA is not a valid value" )

    # validate value type and range
    fnCurr <- dfValid$sFnValidate[rr]
    if( is.na(fnCurr)) next
    minCurr <- dfValid$min[rr]
    maxCurr <- dfValid$max[rr]
    vValue <- do.call( fnCurr, list(dfPed[,sCol],minCurr,maxCurr) )
    if( is.character(vValue) ) stop( "Problem with pedigree column '",sCol,"': ", vValue )

    # update ped (can change type)
    dfPed[,sCol] <- vValue
  }

  # validate - famid
  sCol <- "famid"
  viuFamId <- unique(dfPed[,sCol])
  if( length(viuFamId)>1 ) stop( "Problem with pedigree column '",sCol,"': ", "Contains more than one value. Only single pedigrees are currently handled")

  # validate - unique family-person ids
  tblId <- table(paste0("famid=",dfPed$famid,", id=",dfPed$id))
  vbxIdsDuplicated <- tblId !=1
  if( sum(vbxIdsDuplicated)){
    sErrMsg <- paste( "Ped file contains duplicate person id(s) -", paste( names(tblId)[vbxIdsDuplicated], collapse=", ") )
    stop(sErrMsg)
  }

  # check all parents have a row in the ped file
  # this is checked for by pedigree() but reported as bad momid - misleading as momid not a ped file column
  # strangely familycheck() doesn't report it as a problem, just reported as an unrelated individual via familycheck$unrelated, again misleading
  sParent <- 'motherid'
  dfMum <- unique(dfPed[dfPed[,sParent]!=0, c('famid',sParent)])
  colnames(dfMum)[colnames(dfMum)==sParent] <- 'id'
  sParent <- 'fatherid'
  dfDad <- unique(dfPed[dfPed[,sParent]!=0, c('famid',sParent)])
  colnames(dfDad)[colnames(dfDad)==sParent] <- 'id'
  dfParent <- rbind( dfMum, dfDad)
  #fnStr(dfParent)
  for(rr in 1:nrow(dfParent)) {
    vbx <- T; for(sCol in colnames(dfParent)) vbx <- vbx & dfPed[,sCol]==dfParent[rr,sCol]
    if(sum(vbx)==1) next

    if(sum(vbx)>1) stop("BUG") # validate
    sId <- apply(dfParent[rr,],1,function(x) { sId <- NULL; for(sCol in names(x)) sId <- paste( sId, sCol,'=',x[sCol] ); sId } )
    stop( "No pedigree information specified for person", sId, ", despite their being a parent")
  }

  # validate - familycheck
  dfFamCheck <- with( dfPed, kinship2::familycheck(famid, id, fatherid, motherid))
  #fnStr(dfFamCheck)

  if( dfFamCheck$unrelated!=0 ) {
    sMsg <- paste( "Pedigree contains at least one individual who is not genetically related to any other pedigree member.",
                   "Usually such an individual is the childless spouse of a family member, not in itself a problem.",
                   "However the pedigree diagram cannot include this person as we don't know to whom they are married, so we disallow such pedigrees" )
    stop( sMsg )
  }

  # check parents old enough to have kids
  iMinParentChildAgeGap <- 10
  if( 'age' %in% colnames(dfPed) ){
    for(rr in 1:nrow(dfParent)) {
      vbx <- T; for(sCol in colnames(dfParent)) vbx <- vbx & dfPed[,sCol]==dfParent[rr,sCol]
      if(sum(vbx)!=1) stop( "BUG") # validate
      ix <- which(vbx)
      ageParent <- dfPed$age[ix]
      if(is.na(ageParent)) next

      if(ageParent >= iMinParentChildAgeGap) next
      sParentId <- apply( dfParent[rr,],1,function(x) { sId <- NULL; for(sCol in names(x)) sId <- paste( sId, sCol,'=',x[sCol] ); sId } )
      stop("Isn't this person rather young to be a parent? -", sParentId )
    }
  }

  # check children younger than parents and by a reasonable amount
  if( 'age' %in% colnames(dfPed) ){
    for(rr in 1:nrow(dfParent)) {
      vbx <- T; for(sCol in colnames(dfParent)) vbx <- vbx & dfPed[,sCol]==dfParent[rr,sCol]
      if(sum(vbx)!=1) stop( "BUG") # validate
      ix <- which(vbx)
      ageParent <- dfPed$age[ix]
      if(is.na(ageParent)) next

      vbxChild <- dfPed$motherid==dfParent$id[rr] | dfPed$fatherid==dfParent$id[rr]
      dfPedChild <- dfPed[vbxChild,]
      vAgeGap <- ageParent - dfPedChild$age
      vbxBadAgeGap <- !(vAgeGap >= iMinParentChildAgeGap | is.na(vAgeGap))
      if( !sum(vbxBadAgeGap) ) next

      sParentId <- apply( dfParent[rr,],1,function(x) { sId <- NULL; for(sCol in names(x)) sId <- paste( sId, sCol,'=',x[sCol] ); sId } )
      sChildId <- apply( dfPedChild[vbxBadAgeGap,colnames(dfParent)][1,],1,function(x) { sId <- NULL; for(sCol in names(x)) sId <- paste( sId, sCol,'=',x[sCol] ); sId } )
      stop( "The ages of the following people seem unreasonable given one is the parent of the other. ", "Parent -",sParentId, ", Child -", sChildId )
    }
  }

  # further validation regarding special relationships, e.g. twins
  oPed <- fnConstructPedigree( dfPed )

  dfPed
}



# don't call library() in a package
suppressPackageStartupMessages( library(kinship2))

#' Constructs a kinship2::pedigree object
#'
#' Constructs a kinship2::pedigree object from ped format style data.frame.
#' Assumes a valid pedigree, i.e. pedigree checked via fnValidatePedigree().
#' Collates twin info.
#' @param dfPed ped format style data.frame
#' @return kinship2:pedigree object
#' @export
fnConstructPedigree <- function( dfPed )
{
  #ZZZZ
  #### special relationships

  # validate MZs same sex
  sHeadingRelation <- "relation.MZ"
  if( sHeadingRelation %in% colnames(dfPed)) {
    viuRelation <- unique(dfPed[,sHeadingRelation])
    for( iRelation in viuRelation ) {
      if( is.na(iRelation) || !iRelation ) next # ignore zero and NA
      vix <- which(dfPed[,sHeadingRelation] == iRelation)
      # validate same sex
      if(length(unique(dfPed$sex[vix]))!=1) {
        sErrMsg <- paste( "Invalid pedigree, MZ twins with different sexes, persons -", paste( dfPed$id[vix], collapse=", ") )
        stop(sErrMsg)
      }
    }
  }

  # validate each set of MZs/DZs are the same age
  if( "age" %in% colnames(dfPed)) {
    vsHeadingRelation <- c("relation.MZ","relation.DZ")
    for( sHeadingRelation in vsHeadingRelation ) {
      if( ! sHeadingRelation %in% colnames(dfPed)) next
      viuRelation <- unique(dfPed[,sHeadingRelation])
      for( iRelation in viuRelation ) {
        if( is.na(iRelation) || !iRelation ) next # ignore zero and NA
        vix <- which(dfPed[,sHeadingRelation] == iRelation)
        # validate same age
        if(length(unique(dfPed$age[vix]))!=1) {
          sErrMsg <- paste( "Invalid pedigree, twins with different ages, persons -", paste( dfPed$id[vix], collapse=", ") )
          stop(sErrMsg)
        }
      }
    }
  }

  # get all MZ pairings
  mRelationMZ <- NULL; mRelationMZ <- fnMakeRelationMatrix( dfPed, "MZ", 1 )
  mRelationMZ
  # get all DZ pairings
  mRelationDZ <- NULL; mRelationDZ <- fnMakeRelationMatrix( dfPed, "DZ", 2 )
  mRelationDZ
  # drop any DZ pairings that are also MZ pairings
  if( !is.null(mRelationMZ) && !is.null(mRelationDZ) ) {
    for(rr in 1:nrow(mRelationMZ)) {
      vbxDrop <- (apply( mRelationDZ[,1:2,drop=F],1, function(x){ all(x==mRelationMZ[rr,1:2]) } ))
      mRelationDZ <- mRelationDZ[!vbxDrop,,drop=F]
    }
  }

  # aggregate special relationships
  mRelationAll <- rbind(mRelationMZ, mRelationDZ)
  #mRelationAll <- mRelationAll[c(1,5,2),]
  #fnStr(mRelationAll)

  # construct pedigree list object
  oPedigreeList <- tryCatch(
    fnConstructPedigreeList( dfPed, mRelationAll )
    ,warning=function(w){ return(w) }
  )
  if( "warning" %in% class(oPedigreeList)) {
    oWarning <- oPedigreeList; oPedigreeList <- NA
    sWrnMsg <- oWarning$message
    if( sWrnMsg=="no non-missing arguments to max; returning -Inf" ) {
      # redo, this time ignoring warning
      oPedigreeList <- suppressWarnings( fnConstructPedigreeList( dfPed, mRelationAll ) )
    } else {
      fnStr(oWarning)
      stop( "In constructing the pedigree object the following warning was raised: ", sWrnMsg )
    }
  }

  # extract pedigree object
  oPedigree <- NULL
  if( class(oPedigreeList) == "pedigree"){
    oPedigree <- oPedigreeList
  } else {
    # get pedigree object of current pedigree
    iFamId <- viuFamId[1]
    oPedigree <- oPedigreeList[1]
  }
  oPedigree
}



# don't call library() in a package
suppressPackageStartupMessages( library(kinship2) )

#### creates a relation matrix from the ped data.frame
# this is required in creation of special relationships (twins) in the pedigree object
fnMakeRelationMatrix <- function( dfPed, sRelationType, iRelationType )
{
  sHeadingRelation <- paste0("relation.",sRelationType)

  if( ! sHeadingRelation %in% colnames(dfPed)) return()

  # collate relationship pairs
  dfRelation <- NULL
  viuRelation <- unique(dfPed[,sHeadingRelation])
  for( iRelation in viuRelation ) {
    if( is.na(iRelation) || !iRelation ) next # ignore zero and NA
    vix <- which(dfPed[,sHeadingRelation] == iRelation)
    vId <- dfPed$id[vix]

    # valid relations come in pairs or more
    if( length(vix)<2 ) stop( paste("Lone twin found, person id -", vId[1]))

    for( r1 in 1:(length(vId)-1) ) for( r2 in (r1+1):length(vId) ) {
      dfRelation <- rbind( dfRelation, data.frame( id1=vId[r1], id2=vId[r2]) )
    }
  }
  if(is.null(dfRelation)) return()

  dfRelation$code = iRelationType
  mRelation <- as.matrix(dfRelation)
  mRelation
}



# constructor of oPedigreeList requires a 4 way branch because of lack of default value handling
fnConstructPedigreeList <- function( dfPed, mRelationAll )
{
  oPedigreeList <- NA
  if( !is.null(mRelationAll)) {
    if( "deceased" %in% colnames(dfPed) ) {
      oPedigreeList <- kinship2::pedigree( id=dfPed$id, dadid=dfPed$fatherid, momid=dfPed$motherid, sex=dfPed$sex, affected=dfPed$affected
                                 , status=dfPed$deceased,
                                 , relation=mRelationAll )
    } else {
      oPedigreeList <- kinship2::pedigree( id=dfPed$id, dadid=dfPed$fatherid, momid=dfPed$motherid, sex=dfPed$sex, affected=dfPed$affected
                                 , relation=mRelationAll )
    }
  } else {
    if( "deceased" %in% colnames(dfPed) ) {
      oPedigreeList <-
        kinship2::pedigree( id=dfPed$id, dadid=dfPed$fatherid, momid=dfPed$motherid, sex=dfPed$sex, affected=dfPed$affected, status=dfPed$deceased )
    } else {
      oPedigreeList <-
        kinship2::pedigree( id=dfPed$id, dadid=dfPed$fatherid, momid=dfPed$motherid, sex=dfPed$sex, affected=dfPed$affected )
    }
  }
  oPedigreeList
}



fnValidateIntegerVector <- function( vValue, iMin=-Inf, iMax=Inf )
{
  if( !all(is.na(vValue)) & !is.numeric(vValue)) return( "Contains non-numeric values")
  viValue <- as.integer(vValue)
  if(any(vValue != viValue & !is.na(viValue))) return( "Non integers are invalid values")
  if(is.na(iMin)) iMin <- -Inf
  if(is.na(iMax)) iMax <- Inf
  if( any((iMin > viValue | iMax < viValue) & !is.na(viValue)) ) return( "Contains values outside the valid range")
  viValue
}
fnValidateNumericVector <- function( vValue, iMin=-Inf, iMax=Inf )
{
  if( !all(is.na(vValue)) & !is.numeric(vValue)) return( "Contains non-numeric values")
  viValue <- as.numeric(vValue)
  if(is.na(iMin)) iMin <- -Inf
  if(is.na(iMax)) iMax <- Inf
  if( any((iMin > viValue | iMax < viValue) & !is.na(viValue)) ) return( "Contains values outside the valid range")
  viValue
}



