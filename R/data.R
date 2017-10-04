#' Definition of Format of Disease Model File
#'
#' Definition of Format of Disease Model File.
#' Essentially a list of parameters, whether they are required or optional, and what values they can take.
#'
#' @format A data frame with 4 variables:
#' \describe{
#'   \item{category}{XXXX}
#'   \item{parameter}{parameter name}
#'   \item{required}{whether the parameter is required or optional}
#'   \item{type}{parameter's value type, e.g. character, integer, proportion, positive, characterVector, table}
#'   ...
#' }
"dfDialectDisFile"


#' 3 Generation Pedigree.
#'
#' 3 Generation Pedigree.
#'
#' @format A data frame with variables corresponding to ped file format:
#' \describe{
#'   \item{famid}{family id}
#'   \item{id}{person id}
#'   \item{fatherid}{father id}
#'   \item{motherid}{mother id}
#'   \item{sex}{sex}
#'   \item{affected}{affected status, 1=affected, 0=unaffected, Na=unknown}
#'   \item{age}{age}
#'   \item{deceased}{1=deceased,0=living}
#'   ...
#' }
"dfPed3Gen1Aff"


#' Disease Model for Major Depression
#'
#' Disease Model for Major Depression.
#'
#' @format A named list containing epidemiological findings defining the model:
#' \describe{
#'   \item{title}{title of the model}
#'   \item{description}{description of the model}
#'   \item{other}{etc.}
#'   ...
#' }
"lDisDepression"


#' Disease Model for Schizophrenia
#'
#' Disease Model for Schizophrenia.
#'
#' @format A named list containing epidemiological findings defining the model:
#' \describe{
#'   \item{title}{title of the model}
#'   \item{description}{description of the model}
#'   \item{other}{etc.}
#'   ...
#' }
"lDisSchizophrenia"
