#' Disease Risk Predictor
#'
#' This package facilitates
#' * developing disease models from epidemiolgical findings
#' * predicting disease risk for the individuals of a family (pedigree), based on these disease models.
#' Only multifactorial dichotomous disorders are currently modelled.
#'
#' @author Desmond Campbell
#'
#' @useDynLib diseaseRiskPredictor, .registration = TRUE
#' @importFrom graphics par plot abline
#' @importFrom Rcpp evalCpp
#' @importFrom graphics frame lines points polygon segments strheight strwidth text


## #name diseaseRiskPredictor
## @docType package



.onAttach <-function(lib,pkg){
  ver <- read.dcf( file.path(lib, pkg, "DESCRIPTION"), "Version")
  title <- read.dcf( file.path(lib, pkg, "DESCRIPTION"), "Title")
  desc <- read.dcf( file.path(lib, pkg, "DESCRIPTION"), "Description")
  sMsg <- paste( title, ver, "loaded\n", desc )
  packageStartupMessage( sMsg )
}
