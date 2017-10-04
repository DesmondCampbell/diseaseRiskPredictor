fnReadNamedValue <- function( sName, vsLines )
{
  #XXXX_SOFT_TAB_PROBLEM sRegExp <- paste0("^",sName,"\\t")
  sRegExp <- paste0("^",sName,"( |\\t)")
  iLine <- grep( sRegExp, vsLines, ignore.case=T )
  if(!length(iLine)) return() # name not found

  if( length(iLine)>1 ) stop( paste0("Multiple name value pairs for - '",sName,"'") )

  # read value from line
  #XXXX_SOFT_TAB_PROBLEM vsParsed <- strsplit(vsLines[iLine], split="\t")[[1]]
  vsParsed <- strsplit(vsLines[iLine], split=c("( |\t)+"))[[1]]
  if(length(vsParsed)==2) return(vsParsed[2]) else
    if(length(vsParsed)==1) return("") else
      stop("Unexpected number of items in parsed line")
}

fnReadHeading <- function( sHeading, vsLines )
{
  sRegExp <- paste0("^",sHeading,"$")
  iLine <- grep( sRegExp, vsLines, ignore.case=T )
  if(!length(iLine)) return() # heading not found
  iLine
}

fnGetLinesAfterHeading <- function( sHeading, vsLines )
{
  iLine <- fnReadHeading( sHeading, vsLines )
  if(is.null(iLine)) return() # heading not found
  if( length(iLine)>1 ) stop( paste0("Multiple headings - '",sHeading,"'") )

  # calc line number of first line of block
  iLineStart <- iLine + 1
  nofLines <- length(vsLines)
  if(iLineStart > nofLines) stop( paste0("No block following heading - '",sHeading,"'") )

  # calc line number of last line of block
  vixBlankLines <- which(grepl( "^( |\\t)*$", vsLines ))
  vbx <- iLine < vixBlankLines
  iLineEnd <- (vixBlankLines[vbx]-1)[1]
  if(is.na(iLineEnd)) iLineEnd <- nofLines

  # return block extent
  c(iLineStart,iLineEnd)
}



fnReadLines <- function( vsLines, vsTokenNameValuePair, vsTokenLinesAfter, vsTokenTableAfter )
{
  #vsLines <- sub( pattern="#.*$", "", vsLines ) # remove comments
  vsLines <- sub( pattern="[\t ]+$", "", vsLines ) # remove trailing spaces

  lOut <- list() # init object that will represent the information gleaned from the file

  for( sToken in vsTokenLinesAfter ){
    viLineStartEnd <- fnGetLinesAfterHeading( sToken, vsLines )
    if( !is.null(viLineStartEnd) ) lOut[[sToken]] <- vsLines[viLineStartEnd[1]:viLineStartEnd[2]]
  }

  for( sToken in vsTokenNameValuePair ){
    value <- fnReadNamedValue( sToken, vsLines )
    if( !is.null(value) ) lOut[[sToken]] <- value
  }

  for( sToken in vsTokenTableAfter ){
    viLineStartEnd <- fnGetLinesAfterHeading( sToken, vsLines )
    if( !is.null(viLineStartEnd) ) {
      txtCon <- textConnection(vsLines[viLineStartEnd[1]:viLineStartEnd[2]])
      dfTbl <- NA; dfTbl <- tryCatch( {
        #XXXX_SOFT_TAB_PROBLEM lOut[[sToken]] <- utils::read.table( file=txtCon, sep="\t", header=T, stringsAsFactors=F, colClasses="character", skip=viLineStartEnd[1]-1, nrows=viLineStartEnd[2]-viLineStartEnd[1] )
        #lOut[[sToken]] <- utils::read.table( file=txtCon, sep="", header=T, stringsAsFactors=F, colClasses="character", skip=viLineStartEnd[1]-1, nrows=viLineStartEnd[2]-viLineStartEnd[1] )
        lOut[[sToken]] <- utils::read.table( file=txtCon, sep="", header=T, stringsAsFactors=F, colClasses="character" )
      }
      , warning=function(w) w$message
      , error=function(e) e$message
      ) # end of tryCatch
      close(txtCon)
      if( is.character(dfTbl)) stop( "Problem with ", sToken, " table : ", dfTbl )
    }
  }

  lOut
}



#' XXXX
#'
#' XXXX
#' @param sFile XXXX
#' @param vsTokenNameValuePair XXXX
#' @param vsTokenLinesAfter XXXX
#' @param vsTokenTableAfter XXXX
#' @return XXXX
#' @export
fnReadFile <- function( sFile, vsTokenNameValuePair, vsTokenLinesAfter, vsTokenTableAfter )
{
  # read file into a vector of lines
  vsLines <- readLines(sFile, n = -1)

  # convert into read object
  lOut <- fnReadLines( vsLines, vsTokenNameValuePair, vsTokenLinesAfter, vsTokenTableAfter )

  lOut
}



fnTypeConvertAndValidateVector <- function( vParamValue, sType )
{
  switch( sType
          , character = {
            vParamValue <- suppressWarnings( as.character(vParamValue) )
            if( any( is.na(vParamValue) | is.null(vParamValue) ) ) stop( "Not text")
          }
          , characterVector = {
            if( "character" != typeof(vParamValue) ) stop( "Not text")
          }
          , table = {
            if( "data.frame" != class(vParamValue) ) stop( "Not table")
          }
          , numeric = {
            vParamValue <- suppressWarnings( as.numeric(vParamValue) )
            if( any( is.na(vParamValue) | is.null(vParamValue) ) ) stop( "Not a number")
          }
          , proportion = {
            vParamValue <- suppressWarnings( as.numeric(vParamValue) )
            if( any( is.na(vParamValue) | is.null(vParamValue) | vParamValue < 0 | vParamValue > 1 ) ) stop( "Not a number between 0 and 1")
          }
          , positive = {
            vParamValue <- suppressWarnings( as.numeric(vParamValue) )
            if( any( is.na(vParamValue) | is.null(vParamValue) | vParamValue < 0 ) ) stop( "Not a positive number")
          }
          , stop( paste( "Unhandled type encountered - ", sType) )
  )
  vParamValue
}



fnTypeConvertRangeCheck <- function( dfDisDef, lDisFile )
{
  # global parameters
  for( rr in 1:nrow(dfDisDef) ) {
    if( dfDisDef$category[rr] != "global" ) next

    sParam <- dfDisDef$parameter[rr]

    # validate required parameters defined
    if( ! sParam %in% names(lDisFile) ){
      if(dfDisDef$required[rr]) stop( paste0("Required parameter not defined - '", sParam, "'"))
      next # nothing more to check as parameter is absent
    }

    # convert value to correct type and validate
    paramValue <- NULL; paramValue <- tryCatch(
      fnTypeConvertAndValidateVector( lDisFile[[sParam]], dfDisDef$type[rr] )
      ,error=function(e){ return(e) }
    )
    if( "error" %in% class(paramValue) ) {
      oError <- paramValue; paramValue <- NULL
      stop( paste0( "Problem value found for '", sParam, "' parameter : ", oError$message ) )
    }
    lDisFile[[sParam]] <- paramValue
  }

  # table parameters
  for( rr in 1:nrow(dfDisDef) ) {
    if( dfDisDef$category[rr] == "global" ) next

    sParam <- dfDisDef$parameter[rr]
    sTable <- dfDisDef$category[rr]

    # if current table not specified in Dis file then ignore this parameter check
    if( !sTable %in% names(lDisFile) ) next

    # validate required parameters defined
    if( ! sParam %in% colnames(lDisFile[[sTable]]) ){
      if(dfDisDef$required[rr]) stop( paste0("Required parameter not defined - table '", sTable, "', column'", sParam, "'"))
      next # nothing more to check as parameter is absent
    }

    # convert value to correct type and validate
    vParamValue <- NULL; vParamValue <- tryCatch(
      fnTypeConvertAndValidateVector( lDisFile[[sTable]][,sParam], dfDisDef$type[rr] )
      ,error=function(e){ return(e) }
    )
    if( "error" %in% class(vParamValue) ) {
      oError <- vParamValue; vParamValue <- NULL
      stop( paste0( "Problem in '", sTable, "' table, '", sParam, "' column : ", oError$message) )
    }
    lDisFile[[sTable]][,sParam] <- vParamValue
  }
  lDisFile
}
