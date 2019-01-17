#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param images An n by V  matrix of voxels inside the brain.
#' @param coefficients An n by V matrix of parameter estimates to subtract from the bootstrapped values.
#' If computing p-values this should be the observed parameter estimates. If computing confidence intervals
#' then this should be 0 (the default).
#' @param X Design matrix for the full model.
#' @param Xred Design matrix for the reduced model.
#' @param W A vector of weights for the regression. Voxel specific weights not accepted
#' @param statistic A function of the form statistic(statimg, mask) that returns a statistic of interest.
#' @keywords statistical parametric map, nonparametric bootstrap
#' @param ... Arguments passed to statistic function.
#' @return Returns a vector of the resulting estimate
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @export
bootStats = function(images,                            # M(n, i)
                     coefficients=0,                    # Maybe M(k-r, i)
                     X,                                 # M(n, k)
                     Xred,                              # M(n, r)
                     statistic=function(stat, ...) stat,# V(i) -> argv -> V(i)
                     ...)                               # argv
{
  rdf    <- nrow(X) - ncol(X)                           # UInt    // n - k

  QR     <- qr(X)                                       # M(n, k)

  # compute numerator of inverse covariance of beta hat
  peind  <- which(!colnames(X) %in% colnames(Xred))     # V(k-r)
  bcoefs <- qr.coef(QR, images)[peind,] - coefficients  # M(k-r, i)
  varX1  <- qr.resid(qr(Xred), X[,peind])               # M(n, k-r)
  varX1  <- t(varX1) %*% varX1                          # M(k-r, k-r)
  stat   <- colSums(bcoefs * (varX1 %*% bcoefs))        # V(i)

  # compute denominator of inverse covariance of beta hat
  images <- qr.resid(QR, images)                        # M(n, i)
  images <- colSums(images^2)/rdf                       # V(i)

  # This is the statistic
  statistic(stat/images, ...)                           # V(i)
}
