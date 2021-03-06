% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstimateRiskSrc.R
\name{fnPredictRisk}
\alias{fnPredictRisk}
\title{Predicts risk in disease pedigree}
\usage{
fnPredictRisk(lPedDis, nofBurnIn = 100, nofDraws = 30000, nofYears = 5,
  bVerbose = F)
}
\arguments{
\item{lPedDis}{disease model and pedigree info}

\item{nofBurnIn}{number of initial iterations of Gibbs sampler to be ignored}

\item{nofDraws}{number of draws from posterior pedigree liability distribution}

\item{nofYears}{number of years for n year risk calculation}

\item{bVerbose}{boolean for reporting plots}
}
\value{
list containing liability draws, risk estimates, etc
}
\description{
A Gibbs sampler is used to draw samples from the pedigree's posterior liability distribution.
From this set of draws risk and n year risk is predicted for each pedigree member.
}
