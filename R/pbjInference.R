#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param statistic A user specified function that takes a RNifti image object and computes a particular statistic of interest.
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
pbjInference = function(statMap, statistic = function(image) max(c(image)), nboot=5000, rboot=stats::rnorm, method=c('nonparametric', 't', 'conditional', 'permutation'), ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')

  sqrtSigma <- statMap$sqrtSigma
  sqrtSigma <- if(is.character(sqrtSigma)) readRDS(sqrtSigma) else sqrtSigma
  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  ndims = length(dim(mask))
  stat = rawstat = stat.statMap(statMap)
  template = statMap$template
  df = sqrtSigma$df
  rdf = sqrtSigma$rdf
  robust = sqrtSigma$robust
  HC3 = sqrtSigma$HC3
  transform = sqrtSigma$transform
  method = tolower(method[1])
  obsstat = statistic(stat, ...)

  dims = dim(sqrtSigma$res)
  n = dims[1]
  V=dims[2]

  boots = list()
  bootdim = dim(rboot(n))


  # If sqrtSigma can be stored and accessed efficiently on disk this can be efficiently run in parallel
  pb = txtProgressBar(style=3, title='Generating null distribution')
  tmp = mask
  if(nboot>0){
  for(i in 1:nboot){
    statimg = pbjBoot(sqrtSigma, rboot, bootdim, robust=robust, method = method, HC3=HC3, transform=transform)
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
