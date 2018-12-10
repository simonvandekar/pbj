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
#' to compute an image summary statistic for each bootstrap.
#' @param ... Arguments passed to statistic function.
#' @return Returns a list of length length(cfts)+2. The first two elements contain
#' statMap$stat and statMap$template. The remaining elements are lists containing the following:
#' \item{pvalues}{A vector of p-values corresponding to the cluster labels in clustermaps.}
#' \item{clustermap}{A niftiImage object with the cluster labels.}
#' \item{pmap}{A nifti object with each cluster assigned the negative log10 of its cluster extent FWE adjusted p-value.}
#' \item{CDF}{A bootstrap CDF.}
#' @export
#' @importFrom RNifti readNifti
#' @importFrom abind abind
npbj = function(images, form, formred, mask, data=NULL, W=NULL, template=NULL, nboot=1000, statistic=sei, ...){
  X = getDesign(form, data=data)
  Xred = getDesign(formred, data=data)
  n = nrow(X)
  if(is.null(W)) W = rep(1, n)
  W = sqrt(W)
  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)

  # load in images
  if(class(images)[1] != 'niftiImage'){
    n=length(images)
    images = as.character(images)
    if(nrow(X)!=n)
      stop('length(images) and nrow(X) must be the same.')
    res = do.call(abind::abind, list(RNifti::readNifti(images), along=4))
    } else {
      n = nrow(X)
      res = images
      rm(images)
    }
    dims = dim(res)

  # load mask
  if(class(mask)[1] !='niftiImage'){
    maskimg=as.character(mask)
    mask = RNifti::readNifti(maskimg)
  }

  # check that first input image and mask dimensions are the same
  ndims = length(dim(mask))
  if(any(dims[1:ndims] != dim(mask))  ){
    stop('images and mask dimensions do not match.\n')
  }

  # check that template and mask dimensions are the same
  if( !is.null(template)){
    if(class(template)[1]!='niftiImage'){
      temp = readNifti(as.character(template))
    } else {
      temp = template
    }
      dims = dim(temp)
      rm(temp)
      if(any(dims[1:ndims] != dim(mask))  ){
        stop('template image and mask dimensions (or pixel dimensions) do not match.\n')
      }
  }

  # subset images to mask
  res = t(apply(res, 4, function(x) x[mask!=0]))

  # Compute coefficients
  QR = qr(X * W)
  coefficients = qr.coef(QR, res * W)[peind,]

  # Perform bootstrap using loop or apply or whatever
  samp = sample(1:nrow(pain$data), replace=TRUE)
  bootStats(images=res[samp,], coefficients=coefficients, mask=mask, X=X[samp,], Xred=Xred[samp,], W=W[samp], statistic=statistic, ...)
}
