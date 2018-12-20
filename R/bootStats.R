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
#' @return Returns a list with the fo:
#' \describe{
#'   \item{stat}{The statistical nifti object. If ncol(X) = ncol(Xred)+1, then this is a Z-statistic map, otherwise it is a chi^2-statistic map.}
#'   \item{sqrtSigma}{The 4d covariance object. This is a V by n matrix R, such that R \%*\% t(R) = Sigma.}
#'   \item{mask}{The input mask.}
#'   \item{template}{The background template used for visualization.}
#'   \item{formulas}{A list containing the full and reduced models.}
#'   \item{robust}{A logical indicating the input setting.}
#'   \item{df}{The numerator degrees of freedom of the test.}
#'   \item{rdf}{The numerator degrees of freedom of the test.}
#' }
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @export
bootStats = function(images, coefficients=0, X, Xred, W=NULL, statistic=function(stat) stat, ...){

  n = nrow(images)
  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)
  if(is.null(W)) W = rep(1, n)

  W = sqrt(W)
  QR = qr(X * W)
  # fit model to all image data
  # p X V
  bcoefs = qr.coef(QR, images * W)[peind,] - coefficients
  # compute the part of the inverse covariance of beta hat
  varX1 = qr.resid(qr(Xred * W), X[,peind] * W)
  varX1 = t(varX1) %*% varX1
  images = qr.resid(QR, images * W)
  # overwrite images with the other part of the inverse covariance of beta hat
  images = colSums(images^2)/rdf
  # diag( t(bcoefs) %*% varX1 %*% bcoefs)
  stat = colSums(bcoefs * (varX1 %*% bcoefs))
  # This is a chi-square statistic
  stat = stat/images # * rdf/df
  # convert to chisquared
  #stat = qchisq(pf(stat, df1=df, df2=rdf), df=df)

  out = statistic(stat, ...)
  return(out)
}
