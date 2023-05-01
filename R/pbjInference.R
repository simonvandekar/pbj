#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from lmPBJ function.
#' @param statistic A user specified function that takes an RNifti image object and computes a particular statistic of interest. There are several provided in the pbj package. See referenced functions below.
#' @param null Perform bootstrap under the null or alternative hypotheses?
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. Should return an n vector. Defaults to Rademacher random variable.
#' @param method Character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#' @param runMode character, that controls output. cdf returns the empirical CDFs, bootstrap returns the bootstrapped statistics as a list.
#' @param progress Logical indicating whether to track progress with a progress bar.
#' @param rdata_rds A character string specifying a .rdata or .rds file to write the output. If specified, resampling is run as background process.
#' @param cft_s A vector of robust effect size index (RESI) cluster forming thresholds
#' @param cft_p A vector of p-value cluster forming thresholds
#' @param cft_chisq A vector of p-value cluster forming thresholds
#' @param ... arguments passed to statistic function.
#'
#' @return Returns the statMap object, with a pbj object added. If runMode=='cdf', the first element is the observed statistic value, and the subsequent elements are the CDFs and ROIs, used for computing adjusted p-values and plotting. If runMode=='bootstrap', the first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @seealso [mmeStat()], [maxima()], and [cluster()] for statistic functions. See [lmPBJ()] to create statMap objects. See [image.statMap()], and [table.statMap()] for producing summaries of the results.
#' @export
#' @details This function runs resampling-based methods to perform inference on topological
#' features of neuroimaging statistics. The topological feature is determined by the
#' `statistic` function. Several exist in the pbj package (see below). The statistic
#' function should take an image as input and compute some topological feature from the
#' image. returned as a list. Multiple topological features can be returned, as in
#' [mmeStat()] and [cluster()]. #' To use default methods the statistic must have a logical
#' `rois` argument that outputs an integer valued image identifying where each topological
#' features is located.
#' @example inst/examples/pbjInference.R
pbjInference = function(statMap, statistic = mmeStat, null = TRUE, nboot=5000, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, method=c('wild', 'permutation', 'nonparametric'), runMode=c('cdf','bootstrap'), progress=FALSE, rdata_rds=NULL, cft_s=NULL, cft_p=NULL, cft_chisq=NULL, ...){
  argsList <- c(as.list(environment()), list(...))
  # CFT passed as p value or effect size converted to chi-squared threshold
  argsList$cft = cft_chisq
  if(!is.null(cft_s)){
   argsList$cft = cft_s^2 * statMap$sqrtSigma$n + statMap$sqrtSigma$df
  } else if(!is.null(cft_p)){
    argsList$cft = qchisq(cft_p, df=statMap$sqrtSigma$df, lower.tail=FALSE)
  }
  argsList = argsList[grep('cft_s|cft_p|cft_chisq', names(argsList), invert=TRUE)]

  if(is.null(rdata_rds)){
    argsList = argsList[grep('rdata_rds', names(argsList), invert=TRUE)]
   do.call(pbj::pbjInferenceFG, argsList)
  } else {
   pbj::pbjInferenceBG(argsList)
  }
}

#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from lmPBJ function.
#' @param statistic A user specified function that takes an RNifti image object and computes a particular statistic of interest. There are several provided in the pbj package. See referenced functions below.
#' @param null Perform bootstrap under the null or alternative hypotheses?
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
#' @importFrom methods formalArgs
#' @seealso [mmeStat()], [maxima()], and [cluster()] for statistic functions. See [lmPBJ()] to create statMap objects. See [image.statMap()], and [table.statMap()] for producing summaries of the results.
#' @export
#' @details This function runs a bootstrap or permutation procedure to perform inference on topological features of neuroimaging statistics.
#' The topological feature is determined by the `statistic` function. Several exist in the pbj package. The statistic function should take an image as input
#' and compute some topological feature from the image. returned as a list. Multiple topological features can be returned, as in [mmeStat()] and [cluster()].
#' To use default methods the statistic must have a logical `rois` argument that outputs an integer valued image identifying where each topological features is located. Details of the resampling procedures are available in Vandekar et al. (2022).
#' @example inst/examples/pbjInference.R
pbjInferenceFG = function(statMap, statistic = mmeStat, null=TRUE, nboot=5000, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, method=c('wild', 'permutation', 'nonparametric'), runMode=c('cdf','bootstrap'), progress=FALSE, ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')
  runMode = tolower(runMode[1])
  method = tolower(method[1])

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
  transform = tolower(sqrtSigma$transform[1])
  method = tolower(method[1])
  # arguments passed to statistic function
  statArgs = list(...)
  # add mask file if it's required for statistic function
  if(any(grepl('mask', formalArgs(statistic)))) statArgs$mask = mask
  statArgs$stat = stat
  obsstat = do.call(statistic, statArgs)
  rois  = if('rois' %in% names(formals(statistic))){
    do.call(statistic, c(statArgs, rois=TRUE))
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
        statArgs$stat[ mask!=0] = statimg
        boots[[i]] = do.call(statistic, statArgs)
        setTxtProgressBar(pb, round(i/nboot,2))
      }
      close(pb)
    }
  } else {
    tmp = mask
    if(nboot>0){
      for(i in 1:nboot){
        statimg = pbjBoot(sqrtSigma, rboot, method = method, null=null)
        statArgs$stat[ mask!=0] = statimg
        boots[[i]] = do.call(statistic, statArgs)
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
                 # This suppress warnings is for when the bootstrap has no clusters (length zero)
                 suppressWarnings(globCDF <- lapply(boots, function(boot){
                   if(sum(sapply(boot, length))==0){
                     NA
                   } else {
                     wecdf(sapply(boot, max))
                   }
                 } ) )
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
               # for(ind in 1:length(obsstat)){
               #   newInds = order(obsstat[[ind]], decreasing=TRUE)
               #   obsstat[[ind]][] = obsstat[[ind]][newInds]
               #   if(!is.null(rois)){
               #     rois[[ind]][,,] = match(rois[[ind]][,,], newInds)
               #     rois[[ind]][is.na(rois[[ind]][,,])] = 0
               #   }
               # }
               list(obsStat=obsstat, margCDF=margCDF, globCDF=globCDF, ROIs=rois)},
             bootstrap=list(obsStat=obsstat, boots=boots) )
  class(out) = c('pbj', 'list')
  statMap$pbj = out
  return(statMap)
}

#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#'
#' @param argsList A list of arguments obtained from pbjInference including `rdata_rds` that specifies a file to save the output.
#' @importFrom callr r_bg
#' @seealso [pbjInference()], [pbjInferenceFG()]
#' @export
#' @details This function performs pbj resampling-based inference running the analyses in the background and saving the results as and rdata or rds object.
#' If the file name ends in rdata it will save the statMap object and the computing time in an rdata file. Otherwise it will save the statMap object as an rds file.
#' Note, the input to this function is different than [pbjInferenceFG()].
pbjInferenceBG = function(argsList){
  rdata_rds = argsList[[grep('rdata_rds',names(argsList) )]]
  argsList = argsList[grep('rdata_rds',names(argsList), invert=TRUE )]
  rcallRes = r_bg(function(argsList, rdata_rds){
    computeTime = system.time(statMap <- do.call(pbj::pbjInferenceFG, argsList) )
    if(grepl('.rds', tolower(rdata_rds))){
      saveRDS(statMap, file=rdata_rds)
    } else {
      save(statMap, computeTime, file=rdata_rds)
    }
  }, args=list(argsList=argsList, rdata_rds=rdata_rds))
  rcallRes
}
