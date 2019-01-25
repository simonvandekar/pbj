
#include <RcppArmadillo.h>

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


// Couple MACRO helpers for debugging
//#define DEBUG_PBJ(x) Rcout << x << "\n"; // Comment in for Debug ON
#define DEBUG_PBJ(x)                     // Comment in for Debug OFF
#define MTYPE(x,y) DEBUG_PBJ("  " << x << " :: M(" << y.n_rows << ", " << y.n_cols << ")")

// [[Rcpp::export]]
arma::mat calcBootStats(
                         arma::mat  images,                // M(n, i)
                         arma::mat  X,                     // M(n, k)
                         arma::mat  Xred,                  // M(n, r)
                         arma::mat  coefficients,          // Maybe M(k-r, i)
                         arma::uvec peind                  // V(k-r)
                       )
{
  arma::mat Q, R, tmp, bcoefs, varX1, identity(X.n_rows, X.n_rows), denom, stat;
  DEBUG_PBJ("\n## calcBootStats\n")

  identity.eye();

  peind -= 1; // Adjust offset base from 1 to 0

  // QR     <- qr(X)                                      # M(n, k)
  DEBUG_PBJ("QR     <- qr(X)")
  arma::qr(Q,R,X);                                  // Q :: M(n, k), R :: M(k, k)

  // compute denominator of inverse covariance of beta hat
  // denom <- qr.resid(QR, images)                        # M(n, i)
  DEBUG_PBJ("denom <- qr.resid(QR, images)")
  //denom = (identity - Q*Q.t())*images;

  // denom <- colSums(denom^2)/(nrow(X) - ncol(X))        # V(i)                    # Need denom::V(i)
  DEBUG_PBJ("denom <- colSums(denom^2)/(nrow(X) - ncol(X))")
  denom  = sum(square((identity - Q*Q.t())*images), 0)/(X.n_rows - X.n_cols);

  //# compute numerator of inverse covariance of beta hat
  //bcoefs <- qr.coef(QR, images)[peind,] - coefficients  # M(k-r, i)               # Reuse Q, R, reuse tmp for bcoefs
  DEBUG_PBJ("bcoefs <- qr.coef(QR, images)[peind,] - coefficients")
  bcoefs = solve(R, Q.t() * images);
  bcoefs = bcoefs.rows(peind);

  DEBUG_PBJ("  -= coefficients :: M(" << coefficients.n_rows << ", " << coefficients.n_cols << ")")
  if(!coefficients.is_empty()) bcoefs-=coefficients;// M(k-r, i)

  DEBUG_PBJ("varX1  <- qr.resid(qr(Xred), X[,peind])");
  X = X.cols(peind);
  arma::qr(Q,R,Xred);
  identity.resize(Q.n_rows, Q.n_rows);
  identity.eye();
  identity -= Q*Q.t();
  varX1  = identity * X;

  DEBUG_PBJ("varX1  <- t(varX1) %*% varX\n")
  //varX1 = varX1.t() * varX1;

  DEBUG_PBJ("stat   <- colSums(bcoefs * (varX1 %*% bcoefs))")
  stat  = sum(bcoefs % (varX1.t() * varX1 * bcoefs), 0);

  return(stat/denom);
}
