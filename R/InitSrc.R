# An init file with useful start up functions etc

# D Campbell 13.9.17



#' My str() function
#'
#' str() with variable naming and flushing of output to file
#' @param x variable to str()
#' @return NULL
#' @export
fnStr <- function( x ) { cat( "VALUE: ", deparse(substitute(x))," = ", sep=""); utils::str(x); utils::flush.console() }
#fnTra <- function( x ) { cat( "TRACE:", x, "\n" ); utils::flush.console() }

#' My cat() function
#'
#' cat() with prefix and flush
#' @param ... variables to cat()
#' @return NULL
#' @export
fnRep <- function( ... ) { cat( "REPORT:", ..., "\n" ); utils::flush.console() }

#' My fix() function
#'
#' Wrapper for fix()
#' @param dfData table to give to fix()
#' @return NULL
#' @export
fnFix <- function( dfData ) utils::fix( dfData )



#### set stringsAsFactors = FALSE within the package but not outside the package
library(hellno)
# import the namespace to have hellno functions override the base functions
#' @import hellno
#options(stringsAsFactors = FALSE) # hurray!!!



# don't call library() in a package
suppressPackageStartupMessages( library(moments) )

#' Calculates descriptive stats on data.frame columns
#'
#' Calculates descriptive stats on data.frame columns
#' @param dfData data-table to be described
#'
#' @return data.frame containing descriptive stats
#' @export
fnSS <- function( dfData )
{
  dfStats <- data.frame( stringsAsFactors=F,
                         variableId= colnames(dfData),
                         mean=apply( dfData, 2, mean, na.rm=T ),
                         var=apply( dfData, 2, stats::var, na.rm=T ), sd=apply( dfData, 2, stats::sd, na.rm=T ),
                         skewness=apply( dfData, 2, moments::skewness, na.rm=T ),
                         min=tryCatch( apply( dfData, 2, min, na.rm=T ), error=function(e)NA, warning=function(w)NA ),
                         max=tryCatch( apply( dfData, 2, max, na.rm=T ), error=function(e)NA, warning=function(w)NA ),
                         count=nrow( dfData ),
                         count.obs=apply( dfData, 2, function(x){sum(!is.na(x))} ),
                         count.na=apply( dfData, 2, function(x){sum(is.na(x))} )
  )
  dfStats
}
