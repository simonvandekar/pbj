#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from lmPBJ function.
#' @param statistic A user specified function that takes an RNifti image object and computes a particular statistic of interest. There are several provided in the pbj package. See referenced functions below.
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. Should return an n vector. Defaults to Rademacher random variable.
#' @param method Character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#' @param runMode character, that controls output. cdf returns the empirical CDFs, bootstrap returns the bootstrapped statistics as a list.
#' @param progress Logical indicating whether to track progress with a progress bar.
#' @param ... arguments passed to statistic function.
#'
#' @return Returns the statMap object, with a pbj object added. If runMode=='cdf', the first element is the observed statistic value, and the subsequent elements are the CDFs and ROIs, used for computing adjusted p-values and plotting. If runMode=='bootstrap', the first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @seealso [mmeStat()], [maxima()], and [cluster()] for statistic functions. See [lmPBJ()] to create statMap objects. See [image.statMap()], and [table.statMap()] for producing summaries of the results.
#' @export
#' @details This function runs a bootstrap or permutation procedure to perform inference on topological features of neuroimaging statistics.
#' The topological feature is determined by the `statistic` function. Several exist in the pbj package. The statistic function should take an image as input
#' and compute some topological feature from the image. returned as a list. Multiple topological features can be returned, as in [mmeStat()] and [cluster()].
#' To use default methods the statistic must have a logical `rois` argument that outputs an integer valued image identifying where each topological features is located. Details of the resampling procedures are available in Vandekar et al. (2022).
pbjInference = function(statMap, statistic = mmeStat, nboot=5000, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, method=c('wild', 'permutation', 'nonparametric'), runMode=c('cdf','bootstrap'), progress=FALSE, ...){
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


  # If sqrtSigma can be stored and accessed efficiently on disk this can be efficiently run in parallel
  if(progress){
    pb = txtProgressBar(style=3, title='Generating null distribution')
    tmp = mask
    if(nboot>0){
      for(i in 1:nboot){
        statimg = pbjBoot(sqrtSigma, rboot, method = method)
        tmp[ mask!=0] = statimg
        boots[[i]] = statistic(tmp, ...)
        setTxtProgressBar(pb, round(i/nboot,2))
      }
      close(pb)
    }
  } else {
    tmp = mask
    if(nboot>0){
      for(i in 1:nboot){
        statimg = pbjBoot(sqrtSigma, rboot, method = method)
        tmp[ mask!=0] = statimg
        boots[[i]] = statistic(tmp, ...)
      }
    }
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
                   if(sum(sapply(boot, length))==0){
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
                 obsstat[[ind]][] = obsstat[[ind]][newInds]
                 rois[[ind]][,,] = match(rois[[ind]][,,], newInds)
                 rois[[ind]][is.na(rois[[ind]][,,])] = 0
               }
               list(obsStat=obsstat, margCDF=margCDF, globCDF=globCDF, ROIs=rois)},
             bootstrap=list(obsStat=obsstat, boots=boots) )
  class(out) = c('pbj', 'list')
  statMap$pbj = out
  return(statMap)
}

#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' This is a wrapper to pbjInference that runs it in the background using the `mmeStat` function for the `statistic` argument.
#'
#' @param statMap statMap object as obtained from lmPBJ function.
#' @param mask The mask file used for the analysis passed to `mmeStat`
#' @param cft The cluster forming threshold to use if CEI or CMI are specified.
#' @param rdata_rds The location the save the results as an rdata file or an rds file as a character.
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. Should return an n vector. Defaults to Rademacher random variable.
#' @param method Character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#' @param runMode character, that controls output. cdf returns the empirical CDFs, bootstrap returns the bootstrapped statistics as a list.
#' @param max Logical indicating whether to do inference on local maxima.
#' @param CMI Logical indicating whether to do inference for cluster masses.
#' @param CEI Logical indicating whether to do inference for cluster extents.
#'
#' @return Returns the statMap object, with a pbj object added. If runMode=='cdf', the first element is the observed statistic value, and the subsequent elements are the CDFs and ROIs, used for computing adjusted p-values and plotting. If runMode=='bootstrap', the first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @seealso [mmeStat()], [maxima()], and [cluster()] for statistic functions. See [lmPBJ()] to create statMap objects. See [image.statMap()], and [table.statMap()] for producing summaries of the results.
#' @export
#' @details
#' For `rdata_rds`, if the string has an RDS extension the `statMap` object is saved in the RDS file. If it has any other extension the `statMap` object and computing time are saved as an rdata file.
pbjInferenceBG = function(statMap, mask, cft, rdata_rds, nboot=5000, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, method='wild', max=FALSE, CMI=FALSE, CEI=TRUE){
 rcallRes = r_bg(function(statMap, nboot, rboot, method, mask, cft, progress, max, CMI, CEI, rdata){
  computeTime = system.time(statMap <- pbj::pbjInference(statMap, nboot = nboot, method=method, mask = mask, cft = cft, max=max, CMI=CMI, CEI=CEI, runMode='cdf'))
  if(grepl('.rds', rdata_rds)){
    saveRDS(statMap, file=rdata_rds)
  } else {
    save(statMap, computeTime, file=rdata_rds)
  }
}, args=list(statMap=statMap, nboot=nboot, rboot=rboot, method=method, mask=mask, cft=cft, max=max, CMI=CMI, CEI=CEI, rdata=rdata))
 rcallRes
}
