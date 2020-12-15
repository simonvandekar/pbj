#' Performs nonParametric Bootstrap Joint (nPBJ) Inference
#'
#' @param images A character vector of subjects' images or a 4D niftiImage object.
#' @param form A formula describing the full model.
#' @param formred A formula describing the reduced model.
#' @param mask A character string of niftiImage object identifying voxels to be included in the analysis.
#' @param data A data frame with the variables in form for associated with subject images.
#' @param W A weight vector for the weighted least squares regression.
#' @param template A template image that will be included in the returned statMap object.
#' @param nboot Number of bootstraps to perform.
#' @param statistic A function of the form statistic(stat, ...) that takes a vector stat and other arguments
#' to compute an image summary statistic for each bootstrap. E.g. \code{\link{npbj.sei}}.
#' @param ... Arguments passed to statistic function.
#' @return Returns a list of length length(cfts)+2. The first two elements contain
#' statMap$stat and statMap$template. The remaining elements are lists containing the following:
#' \item{pvalues}{A vector of p-values corresponding to the cluster labels in clustermaps.}
#' \item{clustermap}{A niftiImage object with the cluster labels.}
#' \item{pmap}{A nifti object with each cluster assigned the negative log10 of its cluster extent FWE adjusted p-value.}
#' \item{CDF}{A bootstrap CDF.}
#' @seealso \code{\link{npbj.sei}}
#' @export
#' @importFrom RNifti readNifti
npbj = function(images, form, formred, mask, data=NULL, W=NULL, template=NULL, nboot=1000, statistic=npbj.sei, ...){


  X = getDesign(form, data=data)
  Xred = getDesign(formred, data=data)
  statmap = lmPBJ(images=images, form=X, formred=Xred, mask=mask, data=data, W=W, template=template, robust=FALSE, transform=FALSE)

  if(class(images)[1] != 'niftiImage'){
    n=length(images)
    images = as.character(images)
    if(nrow(X)!=n)
      stop('length(images) and nrow(X) must be the same.')
    res = simplify2array(RNifti::readNifti(images))
  } else {
    n = nrow(X)
    res = images
    rm(images)
  }

  # load mask
  if(class(mask)[1] !='niftiImage'){
    maskimg=as.character(mask)
    mask = RNifti::readNifti(maskimg)
  }

  # load images
  res = t(apply(res, 4, function(x) x[mask!=0]))

  cat(format(Sys.time(), "%a %b %d %X %Y"), ' Running bootstrap\n')

  # Apply weighting before bootstrap
  if(!is.null(W))
  {
    W      <- sqrt(W)
    X      <- X * W
    Xred   <- Xred * W
    images <- images * W
  }

  result = lapply(1:nboot, function(ind){
    samp = sample(1:nrow(X), replace=TRUE)
    bootStats(images=res[samp,], X=X[samp,, drop=FALSE], Xred=Xred[samp,,drop=FALSE], coefficients=statmap$coef, statistic=statistic, mask=mask, df=statmap$df, ...)
  }
  )

  cat(format(Sys.time(), "%a %b %d %X %Y"), ' Finished bootstrap\n')

  # output will depend on function "statistic"
  result = do.call(rbind, result)
  result = list(statmap, stats=result)
}
