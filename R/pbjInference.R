#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param statistic A user specified function that takes a RNifti image object and computes a particular statistic of interest.
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. Should return an n vector. Defaults to Rademacher random variable.
#' @param method Character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#' @param runMode character, that controls output. cdf returns the empirical CDFs, bootstrap returns the bootstrapped statistics as a list.
#' @param ... arguments passed to statistic function.
#'
#' @return Returns a list. if runMode=='bootstrap', the first element is the observed statistic value and the second is a list of the boostrap values. If runMode=='cdf', the first element is the observed statistic value, and the subsequent elements are the CDFs and ROIs, used for computing adjusted p-values and plotting.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @export
pbjInference = function(statMap, statistic = function(image) max(c(image)), nboot=5000, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, method=c('wild', 'permutation', 'nonparametric'), runMode=c('bootstrap', 'cdf'), ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')
  runMode = tolower(runMode[1])

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
  rois  = if('rois' %in% names(formals(statistic))){
    statistic(stat, ..., rois=TRUE)
  } else {
    NULL
  }

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

  out=switch(runMode,
             cdf={
               if(is.list(boots[[1]])){
                 # reorders list
                 boots = apply(do.call(rbind, boots), 2, as.list)
                 # add in observed values (plus a constant) to bootstraps. This avoids zero p-values from the empirical CDF
                 #boots = mapply(c, boots, lapply(obsstat, function(x) list(x+1)), SIMPLIFY=FALSE )
                 margCDF = lapply(boots, function(boot){
                   Cs = sapply(boot, length)
                   if(sum(Cs)==0){
                     NA
                   } else {
                     wecdf(unlist(boot), rep(1/Cs, Cs) )
                   }
                 } )
                 globCDF = lapply(boots, function(boot){
                   Cs = sapply(boot, length)
                   if(sum(Cs)==0){
                     NA
                   } else {
                     wecdf(sapply(boot, max))
                   }
                 } )
               } else {
                 Cs = sapply(boots, length)
                 if(sum(Cs)==0){
                   margCDF <- globCDF <- 0
                 } else {
                   margCDF = wecdf(unlist(boots), rep(1/Cs, Cs) )
                   globCDF = wecdf(sapply(boots, max))
                 }
               }

               # reindex ROIs and obsStat
               for(ind in 1:length(obsstat)){
                 newInds = order(obsstat[[ind]], decreasing=TRUE)
                 obsstat[[ind]] = obsstat[[ind]][newInds]
                 rois[[ind]][,,] = match(rois[[ind]][,,], newInds)
                 rois[[ind]][is.na(rois[[ind]][,,])] = 0
               }
               list(obsStat=obsstat, margCDF=margCDF, globCDF=globCDF, ROIs=rois)},
             bootstrap=list(obsStat=obsstat, boots=boots) )
  return(out)
}
