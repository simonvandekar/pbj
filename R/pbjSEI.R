#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Spatial Extent Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param cfts.s Numeric vector of cluster forming thresholds to use based on the robust effect size index.
#' @param cfts.p Numeric vector of clusterforming thresholds to use based on p-value thresholding. Consistent
#' with other imaging software these thresholds are one tailed (e.g. p<0.01 imples z>2.32).
#' @param nboot Number of bootstrap samples to use.
#' @param kernel Kernel to use for computing connected components. box is
#'  default (26 neighbors), but diamond may also be reasonable. argument to mmand::shapeKernel
#' @param rboot Function for generating random variables. See examples.
#' @param method character method to use for bootstrap procedure.
#' @param debug Returns extra output for statistical debugging.
#'
#' @return Returns a list of length length(cfts)+4. The first four elements contain
#' statMap$stat, statMap$template, statMap$mask, and statMap$df. The remaining elements are lists containing the following:
#' \item{pvalues}{A vector of p-values corresponding to the cluster labels in clustermaps.}
#' \item{clustermap}{A niftiImage object with the cluster labels.}
#' \item{pmap}{A nifti object with each cluster assigned the negative log10 of its cluster extent FWE adjusted p-value.}
#' \item{boots}{The bootstrap values.}
#' \item{obs}{The size of the observed contiguous clusters.}
#' @export
#' @importFrom stats ecdf qchisq rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti writeNifti updateNifti
#' @importFrom mmand shapeKernel
pbjSEI = function(statMap, cfts.s=c(0.1, 0.25), cfts.p=NULL, nboot=5000, kernel='box', rboot=stats::rnorm, method=c('robust', 't', 'regular')){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')
  debug = getOption('pbj.debug', default=FALSE)

  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  ndims = length(dim(mask))
  rawstat = stat.statMap(statMap)
  template = statMap$template
  sqrtSigma = statMap$sqrtSigma
  df = statMap$sqrtSigma$df
  rdf = statMap$sqrtSigma$rdf
  n = statMap$sqrtSigma$n
  robust = statMap$robust
  HC3 = statMap$HC3
  transform = statMap$transform
  stat = rawstat
  method = tolower(method[1])

  if(!is.null(cfts.p)){
    es=FALSE
    cfts = cfts.p
    cfts.s = NULL
  } else {
    es=TRUE
    cfts = cfts.s
  }
  cftsnominal = cfts

  if(es){
    ts = cfts^2 * n + df
  } else {
    ts = qchisq(cfts, df, lower.tail=FALSE)
  }
  stat = rawstat
  # ts are chi-squared statistic thresholds

  tmp = mask
  ndims = length(dim(mask))
  tmp = lapply(ts, function(th){ tmp[ mask!=0] = (stat[mask!=0]>th); tmp})
  k = mmand::shapeKernel(ndims, ndims, type=kernel)
  clustmaps = lapply(tmp, function(tm, mask) {out = mmand::components(tm, k); out[is.na(out)] = 0; RNifti::updateNifti(out, mask)}, mask=mask)
  ccomps = lapply(tmp, function(tm) table(c(mmand::components(tm, k))) )

  rm(statMap)


  boots = matrix(NA, nboot, length(cfts))
  if(debug) statmaps = rep(list(NA), nboot)

  pb = txtProgressBar(style=3, title='Generating null distribution')
  for(i in 1:nboot){
    tmp = mask
    boot = rboot(n)
    bootdim = dim(boot)
    statimg = pbjBoot(sqrtSigma, rboot, bootdim, robust=robust, method = method, HC3=HC3, transform=transform)
      if(debug) statmaps[[i]] = statimg
      tmp = lapply(ts, function(th){ tmp[ mask!=0] = (statimg>th); tmp})
      boots[i, ] = sapply(tmp, function(tm) max(c(table(c(mmand::components(tm, k))),0), na.rm=TRUE))
      setTxtProgressBar(pb, round(i/nboot,2))
    }
    close(pb)
  # } else { # Support for Shared memory
  #
  #   cat('Generating null distribution in parallel.\n')
  #
  #   jobs <- lapply(1:nboot, function(i) {
  #     mcparallel({
  #       require(bigalgebra)
  #       tmp     <- mask
  #       S       <- matrix(rnorm(r*df), r, df)
  #       statimg <- rowSums((sqrtSigma %*% S)^2)
  #       tmp     <- lapply(ts, function(th){ tmp[ mask!=0] = (statimg>th); tmp})
  #       sapply(tmp, function(tm) max(c(table(c(mmand::components(tm, k))),0), na.rm=TRUE))
  #       return(123)
  #     }, detached=FALSE)
  #
  #   })
  #   Fs = mccollect(jobs, wait=TRUE)
  #   Fs = mccollect(jobs)
  #   browser()
  #   Fs = do.call(rbind, Fs)
  #   cat('Completed generation of null distribution.\n')
  # }
  rm(sqrtSigma) # Free large big memory matrix object

  # add the stat max
  ccomps = lapply(ccomps, function(x) if(length(x)==0) 0 else x)
  pvals = lapply(1:length(cfts), function(ind) colMeans(outer(boots[,ind], ccomps[[ind]], FUN = '>') ) )
  names(pvals) = paste('cft', ts, sep='')
    pmaps = lapply(1:length(ts), function(ind){ for(ind2 in 1:length(pvals[[ind]])){
      clustmaps[[ind]][ clustmaps[[ind]]==ind2] = -log10(pvals[[ind]][ind2])
    }
      clustmaps[[ind]]
    } )
  names(pvals) <- names(pmaps) <- names(clustmaps) <- if(es) paste0('cft.s', cftsnominal) else paste0('cft.p', cftsnominal)
  out = list(pvalues=pvals, clustermap=clustmaps, pmap=pmaps, boots=apply(boots, 2, list), obs=ccomps)
  # changes indexing order of out
  out = apply(do.call(rbind, out), 2, as.list)
  out = c(stat=list(rawstat), template=list(template), mask=list(mask), df=list(df), out)
  if(debug) out$boots = statmaps
  class(out) = c('pbj', 'list')
  return(out)
}
