#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param images An n by V  matrix of voxels inside the brain.
#' @param coefficients An n by V matrix of parameter estimates.
#' @param mask A nifti image that defines where the columns of images comes from.
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
bootStats = function(images, coefficients, mask, X, Xred, W=NULL, statistic=function(stat) stat, ...){

  n = nrow(images)
  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)

  W = sqrt(W)
  QR = qr(X * W)
  # fit model to all image data
  if(df==1){
    num = qr.coef(QR, images * W)[peind,]
    seX1 = sqrt(chol2inv(qr.R(QR))[peind,peind])
    images = t(qr.resid(QR, images * W))
    rm(QR)
    # overwrite images
    images = rowSums(images^2)
    stat = sqrt(images/rdf)
    stat = num/stat / seX1
    # convert to z-statistics
    stat = qnorm(pt(stat, df=rdf))
  } else {
    # p X V
    bcoefs = qr.coef(QR, images * W)[peind,] - coefficients
    # compute the part of the inverse covariance of beta hat
    varX1 = chol2inv(qr.R(qr( qr.resid(qr(Xred * W), X[,peind] * W)) ) )
    images = t(qr.resid(QR, images * W))
    # overwrite images with the other part of the inverse covariance of beta hat
    images = rowSums(images^2)/rdf
    # diag( t(bcoefs) %*% varX1 %*% bcoefs)
    stat = colSums((varX1 %*% bcoefs) * bcoefs)
    # This is a chi-square statistic
    stat = stat/images # * rdf/df
    # convert to chisquared
    #stat = qchisq(pf(stat, df1=df, df2=rdf), df=df)
  }

  out = statistic(stat, mask, ...)
  return(out)
}
