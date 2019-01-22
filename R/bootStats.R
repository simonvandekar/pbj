#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param images An n by V  matrix of voxels inside the brain.
#' @param X Design matrix for the full model.
#' @param Xred Design matrix for the reduced model.
#' @param W A vector of weights for the regression. Voxel specific weights not accepted
#' @param coefficients An n by V matrix of parameter estimates to subtract from the bootstrapped values.
#' If computing p-values this should be the observed parameter estimates. If computing confidence intervals
#' then this should be NULL (the default).
#' @param statistic A function of the form statistic(statimg, mask) that returns a statistic of interest.
#' @keywords statistical parametric map, nonparametric bootstrap
#' @param ... Arguments passed to statistic function.
#' @return Returns a vector of the resulting estimate
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @export
bootStats = function(images,                            # M(n, i)
                     coefficients=NULL,                 # Maybe M(k-r, i)
                     X,                                 # M(n, k)
                     Xred,                              # M(n, r)
                     statistic=function(stat, ...) stat,# V(i) -> argv -> V(i)
                     ...)                               # argv
{
  peind  <- which(!colnames(X) %in% colnames(Xred))     # V(k-r)                  # Need peind :: V(k-r)
  # Call cpp instead
  return(statistic(calcBootStats(images, X, Xred, coefficients, peind), ...))

  QR     <- qr(X)                                       # M(n, k)                 # Need Q, R

  # compute denominator of inverse covariance of beta hat
  denom <- qr.resid(QR, images)                         # M(n, i)                 # Need tmp :: M(n,i)
  denom <- colSums(denom^2)/(nrow(X) - ncol(X))         # V(i)                    # Need denom::V(i)

  # compute numerator of inverse covariance of beta hat
  bcoefs <- qr.coef(QR, images)[peind,] - coefficients  # M(k-r, i)               # Reuse Q, R, reuse tmp for bcoefs
  varX1  <- qr.resid(qr(Xred), X[,peind])               # M(n, k-r)               # Need varX1 :: M(n, k-r)
  varX1  <- t(varX1) %*% varX1                          # M(k-r, k-r)
  stat   <- colSums(bcoefs * (varX1 %*% bcoefs))        # V(i)                    # May reuse denom unnamed

  # This is the statistic
  statistic(stat/denom, ...)                            # V(i)                    # Total direct alloc, Q, R, M(n,i), V(i), V(k-r), M(n, k-r)
}
