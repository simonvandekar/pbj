# functions for pbj objects

#' @export
print.pbj <- function(x, ...)
{
  summary(x, ...)
}

cat0 <- function(...) cat(..., sep='')

#' @export
#' @include statmap.R
#' @importFrom stats quantile
summary.pbj <- function(object, ...)
{
  cat0(
    "\nContents:\n",
    statInner("  Stat:       ", object$stat),
    statInner("  Template:   ", object$template)
  )
  
  for(cft in names(object)[ ! names(object) %in% c('stat', 'template') ]){
    cat0('\n', cft, ':\n')
    
    cat0("  P-Values:\n")
    print(quantile(object[[cft]]$pvalues))
    
    cat0(statInner("  Cluster Map: ", object[[cft]]$clustermap),
         statInner("  P Map:       ", object[[cft]]$pmap)
    )
  }
}

#' image a pbj object
#' 
#' See image.statMap for additional arguments
#' 
#' @export
#' @param x pbj object to create images for
#' @param alpha numeric; threshold to apply for threshold mask of significant voxels
#' @param ... Arguments passed to image.statMap
#' @importFrom graphics image
#' @importFrom graphics par
image.pbj <- function(x, alpha=0.05, ...)
{
  for(cft in names(x)[ ! names(x) %in% c('stat', 'template') ]){
    stat = x$stat
    x[[cft]]
    # mask stat image with significant voxels
    stat[ abs(x[[cft]]$pmap) < (-log10(alpha) ) ] = 0
    # create a barebones statmap object
    statmap = list(stat=stat, sqrtSigma=NULL, mask=NULL, template=x$template, formulas=NULL, robust=NULL) 
    class(statmap) = "statMap"
    # call image.statMap
    image(statmap, thresh=qnorm(1-as.numeric(gsub("[^0-9\\.]", "", cft)) )  )
  }
}

#' Write a pbj object to disk
#' 
#' Write a pbj object to disk in parts
#' 
#' @param x pbj object to write
#' @param outdir output directory to write pbj pieces
#' @param ... additional arguments; unused.
#' @export
write.pbj <- function(x, outdir, ...)
{
  statimg = file.path(outdir, 'stat.nii.gz')
  if(!file.exists(statimg)){ 
    if(is.character(x$stat)){
      file.copy(x$stat, statimg)
    } else {
      writeNifti(x$stat, statimg)
    }
  }
  for(cft in names(x)[ ! names(x) %in% c('stat', 'template') ]){
    pmapimg = file.path(outdir, paste0('pbj_sei_log10p_', cft, '.nii.gz'))
    clustmapimg = file.path(outdir, paste0('pbj_sei_clust_', cft, '.nii.gz'))
    writeNifti(x[[cft]]$pmap, pmapimg)
    writeNifti(x[[cft]]$clustermap, clustmapimg)
  }
}
