#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>
#include <algorithm>
#include <random>
#include <ctime>
// [[Rcpp::plugins(openmp)]]

// see http://stackoverflow.com/questions/29311481/sourcecpp-upcoming-iso-standard
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// returns cumulative maximums only in the intervals defined by maxes
arma::vec cummaxk(arma::vec x, arma::vec maxes){
  int k = maxes.n_elem;
  arma::vec out(k);
  out(0) = x(0);
  for(int i = 1; i < k; i++){
    // doesn't matter that there's overlap on the edges b/c we will take cummax
    out(i) = max( x.subvec(maxes(i-1)-1, maxes(i)-1) );
    out(i) = (out(i)<out(i-1))?out(i-1):out(i);
  }
  return out;
}

// returns cumulative maximums
arma::vec cummax(arma::vec x){
  arma::vec out(x);
  out(0) = x(0);
  for(int i = 1; i < x.n_elem; i++){
    out(i) = (out(i)<out(i-1))?out(i-1):out(i);
  }
  return out;
}

arma::vec randChi( arma::mat V, int df){
  arma::mat X = arma::randn<arma::mat>(V.n_cols, df); // n x m_1
  arma::vec num(V.n_rows);
  num = sum(square(V * X), 1);
  return(num);
}

// This function assumes that all rows of M are unit vectors
// It does not check that this is true.
// mu is a T or chisquare statistical map
// f is a positive constant that is a function of the minimum R^2 required
// df is numerator dof for the test
// [[Rcpp::export]]
arma::mat pbjES(arma::vec mu, arma::mat M, signed long chsq, int df, int nboot) {
  arma::mat O(nboot, 2);
  for(int i = 0; i < nboot; i++) {
    // The expected value of the F-stat is \beta^T\Sigma^{1}_\beta \beta + df
    //Rcout << randChi( M, df).size() << '\n' << mu.size() << '\n';
    // subtract out df because it gets added in by randChi
    arma::vec U = randChi( M, df) + mu - df;
    O(i, 0) = arma::min(U( find(mu>chsq) ));
    O(i, 1) = arma::max(U( find(mu<=chsq) ));
  }
  return O;
}

// [[Rcpp::export]]
arma::mat pbjESboundary(arma::mat M, int nboot) {
  arma::mat O(nboot, 1);
  for(int i = 0; i < nboot; i++) {
    // The expected value of the F-stat is \beta^T\Sigma^{1}_\beta \beta + df
    //Rcout << randChi( M, df).size() << '\n' << mu.size() << '\n';
    // subtract out df because it gets added in by randChi
    arma::mat X = arma::randn<arma::mat>(M.n_cols, 1);
    arma::vec U = arma::abs(M * X);
    O(i) = arma::max(U);
  }
  return O;
}

