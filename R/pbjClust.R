#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Cluster Extent Inference
#'
#' @param stat observed statistical map e.g. as obtained from computeStats.
#' @param res Character of 4d residuals in nii or nii.gz format, or residual
#'  matrix as obtained from computeStats.
#' @param mask Character mask file location.
#' @param statoutfiles Character basename for nii.gz files to write out
#'  -log10(p) maps. These are thresholded maps, where the value in each cluster
#'  is the cluster adjusted p-value.
#' @param df Degrees of freedom of test statistics. This is the numerator
#'  degrees of freedom of the F-statistic. If you are passing a Z-statistic
#'  image set df=0.
#' @param rdf Residual degrees of freedom. This is the denominator degrees of
#'  freedom of the F-statistic. This is an optional parameter. If rdf<<n then
#'  this might save time, otherwise, it is ok to leave NULL.
#' @param cfts Numeric vector of cluster forming thresholds to use. These are
#'  single-tailed for chi-squared (or F) statistics or two tailed for
#'  Z-statistics probabilities.
#' @param nboot Number of bootstraps to perform.
#' @param kernel Kernel to use for computing connected components. Box is
#'  default (26 neighbors), but Diamond may also be reasonable.
#'
#' @return Returns a list with the following values:
#' \item{pvalues}{A list of length length(cfts) where each element is a vector of p-values corresponding to the cluster labels in clustermaps.}
#' \item{clustermaps}{A list of length length(cfts) where each element is a nifti object with the cluster labels.}
#' \item{pmaps}{A list of length length(cfts) where each element is a nifti object with each cluster assigned its cluster extent FWE adjusted p-value.}
#' @export
# @examples
pbjClust = function(stat=NULL, res=NULL, mask=NULL, statoutfiles=NULL, df=0, rdf=NULL, cfts=c(0.01, 0.005), nboot=5000, kernel='box'){

  if(is.null(mask))
    stop('Must specify mask.')

  if(is.null(stat))
    stop('Must pass stat map.')

  if(is.null(res))
    stop('Must specify res matrix or images.')

  if(is.character(mask))
    mask=readNifti(mask)

  if(is.character(stat))
    stat = readNifti(stat)

  if(df==0){
    ts = qchisq(cfts, 1, lower.tail=FALSE)
    sgnstat = sign(stat)
    stat = stat^2
  } else {
    ts = qchisq(cfts, df, lower.tail=FALSE)
  }

  tmp = mask
  tmp = lapply(ts, function(th){ tmp[ mask==1] = (stat[mask==1]>th); tmp})
  k = mmand::shapeKernel(3, 3, type=kernel)
  clustmaps = lapply(tmp, function(tm) {out = mmand::components(tm, k); out[is.na(out)] = 0; RNifti::updateNifti(out, mask)})
  stat = lapply(tmp, function(tm) table(c(mmand::components(tm, k))) )

  if(is.character(res)){
    res = readNifti(res)
    res = apply(res, 4, function(x) x[mask==1])
  }
  n = ncol(res)
  if(is.null(rdf)) rdf=n

  ssqs = sqrt(rowSums(res^2))
  if( any(ssqs != 1) ){
    res = sweep(res, 1, ssqs, '/')
  }
  p=nrow(res)
  # This SVD is not strictly necessary
  # If rdf<<n than n and there is a large number of simulations it might save time
  r = min(rdf, p)
  if(rdf != n){
    res = svd(res, nu=r, nv=0)
    res = sweep(res$u, 2, res$d[1:r], "*")
  }

  Fs = matrix(NA, nboot, length(cfts))
  pb = txtProgressBar(style=3, title='Generating null distribution')
  for(i in 1:nboot){
    tmp = mask
    S = matrix(rnorm(r*df), r, df)
    statimg = rowSums((res %*% S)^2)
    tmp = lapply(ts, function(th){ tmp[ mask==1] = (statimg>th); tmp})
    Fs[i, ] = sapply(tmp, function(tm) max(c(table(c(mmand::components(tm, k))),0), na.rm=TRUE))
    setTxtProgressBar(pb, round(i/nboot,2))
  }
  close(pb)
  # add the stat max
  stat = lapply(stat, function(x) if(length(x)==0) 0 else x)
  Fs = rbind(Fs, sapply(stat, max)+0.01) # + 0.01 to make it larger than observed data
  # compute empirical CDFs
  Fs = apply(Fs, 2, ecdf)
  pvals = lapply(1:length(cfts), function(ind) 1-Fs[[ind]](stat[[ind]]) )
  names(pvals) = paste('cft', ts, sep='')
  if(df>0){
    pmaps = lapply(1:length(ts), function(ind){ for(ind2 in 1:length(pvals[[ind]])){
      clustmaps[[ind]][ clustmaps[[ind]]==ind2] = -log10(pvals[[ind]][ind2])
    }
      clustmaps[[ind]]
    } )
  } else {
    pmaps = lapply(1:length(ts), function(ind){ for(ind2 in 1:length(pvals[[ind]])){
      clustmaps[[ind]][ clustmaps[[ind]]==ind2] = -log10(pvals[[ind]][ind2]) * sgnstat
    }
      clustmaps[[ind]]
    } )
  }
  names(pvals) <- names(pmaps) <- names(clustmaps) <- paste('cft', cfts, sep='')

  if(!is.null(statoutfiles)){
    for(cft in cfts){
      RNifti::writeNifti(pmaps[[ paste('cft', cft, sep='') ]], file=paste(gsub('.nii.gz', '', statoutfiles), '_cft', cft, '.nii.gz', sep='') )
    }
  }

  list(pvalues=pvals, clustermaps=clustmaps, pmaps=pmaps)
}
