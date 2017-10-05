#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector fnTimesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
fnTimesTwo(42)
*/




#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix na_matrix(const int nrow, const int ncol){
  NumericMatrix m(nrow, ncol) ;
  std::fill( m.begin(), m.end(), NumericVector::get_na() ) ;
  return m;
}

// [[Rcpp::export]]
double scalarN(const double m, const double s) {
  return R::rnorm(m, s);
}

// [[Rcpp::export]]
double scalarU(const double m, const double s) {
  return R::runif(m, s);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cppFnGibbsSampler( Rcpp::DataFrame dfPed, const Rcpp::NumericMatrix mPriorMean, const Rcpp::List lmRegCoeff, const Rcpp::NumericVector vPostStdDev,
                                       const Rcpp::NumericVector vDrawInit,
                                       const int nofBurnIn, const int nofIter, const int nofVariables, const Rcpp::CharacterVector vsVariables,
                                       const Rcpp::CharacterVector vPerson, const Rcpp::LogicalVector vbxLia, const int iMaxAttempts = 500000 ) {

  // specify starting point in parameter space for Gibbs Sampler
  NumericVector vDraw = vDrawInit;
  vDraw.attr("names") = VECTOR_ELT(mPriorMean.attr("dimnames"), 0);
  NumericVector vPriorMean = mPriorMean(_, 0);
  vPriorMean.attr("names") = VECTOR_ELT(mPriorMean.attr("dimnames"), 0);
  const CharacterVector vsPerson = as<CharacterVector>(dfPed["id"]);
  const NumericVector vExprProp = dfPed["expressedProportionOfLifetimeRisk"];
  const IntegerVector vAffected = dfPed["affected"];
  const NumericVector vThr = dfPed["thr"];

  // init draws repository
  Rcpp::NumericMatrix mDraws = na_matrix(nofIter, nofVariables);
  Rcpp::colnames(mDraws) = vsVariables;

  for (int ii = -nofBurnIn; ii < nofIter; ++ii) {
    if ( ii % 1000 == 0 ) {
      Rcout << "INFO: iteration " << ii << " of " << nofIter << "\n";
    }

    for (int vv = 0; vv < vsVariables.length(); ++vv) {

      // the variable to be conditioned
      const char * sB = vsVariables[vv];

      // vsA -- the set of variables to be conditioned on
      LogicalVector vbxA(vsVariables.length());
      for( int qq = 0; qq < vsVariables.length(); ++qq ) vbxA[qq] = strcmp(vsVariables[qq], sB);
      const CharacterVector vsA = vsVariables[vbxA];

      // calc posterior distribution
      // mean
      const NumericVector vRegCoeffA = lmRegCoeff[sB];
      const NumericVector vDrawA = vDraw[vsA];
      const NumericVector vPMA = vPriorMean[vsA];
      double meanB = vPriorMean[sB];
      for( int aa = 0; aa < vsA.length(); ++aa ) meanB += vRegCoeffA[aa] * (vDrawA[aa] - vPMA[aa]);
      // std dev
      const double stdDevB = vPostStdDev[sB];

      if ( !vbxLia[vv] ) {
        // draw sample of gaussian variable
        vDraw[sB] = scalarN(meanB, stdDevB);

      } else {
        // draw sample of total liability conditioned on affection status and expressed proportion of lifetime risk

        // find index in pedigree for current person
        const char * sPerson = vPerson[vv];
        int pp = -99;
        for (int pcounter = 0; pcounter < vsPerson.length(); ++pcounter) {
          if (!strcmp(sPerson, vsPerson[pcounter])) {
            if (pp >= 0) {
              char sErrMsgP[500];
              sprintf( sErrMsgP, "Expected exactly one person in pedigree with id = '%s'.", sPerson );
              stop(sErrMsgP);
            }
            pp = pcounter;
          }
        }

        // draw a sample that matches the person's observed affection status
        double liaB;
        int ss = 0;
        while(TRUE) {
          // draw sample
          liaB = scalarN(meanB, stdDevB);

          // draw affection status based on total liability
          bool bAffB = liaB > vThr[pp];
          if (bAffB) bAffB = scalarU(0, 1) < vExprProp[pp]; // if above threshold draw affection status depending on expressed proportion of lifetime risk

          // check whether the drawn affection status matches observed
          if ( IntegerVector::is_na(vAffected[pp]) ) break; // their affection status is unknown so any liability matches
          if ( vAffected[pp] == 0 && !bAffB ) break;
          if ( vAffected[pp] == 1 &&  bAffB ) break;

          // report problem if maximum number of attempts exceeded
          if( ++ss > iMaxAttempts ) {
            char sErrMsg[500];
            sprintf( sErrMsg, "The pedigree seems very unlikely given the disease model. Was unable to draw a liability sample for person %s despite %d attempts.", sPerson, iMaxAttempts );
            stop(sErrMsg);
          }
        }
        // record liability draw
        vDraw[sB] = liaB;
      }
    }
    if( ii < 0 ) continue; // burn in

    // record draw for pedigree
    for( int qq = 0; qq < vsVariables.length(); ++qq ) mDraws(ii,qq) = vDraw[qq];
  }

  return mDraws;
}
