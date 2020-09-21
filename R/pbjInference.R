#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param statistic A user specified function that takes a RNifti image object and computes a particular statistic of interest.
#' @param randomX logical, bootstrap samples of X as well?
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. See examples.
#' @param method character method to use for bootstrap procedure.
#' @param ... arguments passed to statistic function.
#'
#' @return Returns a list of length 2. The first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @export
pbjInference = function(statMap, statistic = function(image) max(c(image)), randomX=FALSE, nboot=5000, rboot=stats::rnorm, method=c('nonparametric', 't', 'conditional', 'permutation'), ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')


  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  ndims = length(dim(mask))
  rawstat = stat.statMap(statMap)
  template = statMap$template
  df = statMap$sqrtSigma$df
  rdf = statMap$sqrtSigma$rdf
  robust = statMap$robust
  HC3 = statMap$HC3
  transform = statMap$transform
  stat = rawstat
  method = tolower(method[1])
  obsstat = statistic(stat, ...)

  sqrtSigma <- if(is.character(statMap$sqrtSigma)) {
    apply(readNifti(statMap$sqrtSigma), ndims+1, function(x) x[mask!=0])
  } else {
    statMap$sqrtSigma
  }
  rm(statMap)

  dims = dim(sqrtSigma$res)
  n = dims[1]
  V=dims[2]
  if(is.null(rdf)) rdf=n


  #if(robust & method[1]!='robust'){
  #  BsqrtInv = matrix(apply(sqrtSigma, 2, function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
  #  sqrtSigma = simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), matrix(sqrtSigma[,ind,], nrow=df, ncol=n, byrow=TRUE))) )
  #  sqrtSigma = aperm(sqrtSigma, c(2,3,1))
  #}

  #sqrtSigma <- as.big.matrix(sqrtSigma)
  boots = list()
  bootdim = dim(rboot(n))

  # if(.Platform$OS.type=='windows')
  # {

  pb = txtProgressBar(style=3, title='Generating null distribution')
  tmp = mask
  if(nboot>0){
  for(i in 1:nboot){
    statimg = pbjBoot(sqrtSigma, rboot, bootdim, randomX=randomX, robust=robust, method = method, HC3=HC3, transform=transform)
    tmp[ mask!=0] = statimg
    boots[[i]] = statistic(tmp, ...)
    setTxtProgressBar(pb, round(i/nboot,2))
  }
  close(pb)
  }
  rm(sqrtSigma) # Free large big memory matrix object

  # add the stat max
  out = list(obsStat=obsstat, boots=boots)
  return(out)
}
