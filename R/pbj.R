# functions for pbj objects

#' @useDynLib pbj, .registration=TRUE
#' @export
#' @method print pbj
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
    statInner("  Template:   ", object$template),
    statInner("  Mask:   ", object$mask)
  )

  for(cft in names(object)[ ! names(object) %in% c('stat', 'template', 'mask', 'df') ]){
    cat0('\n', cft, ':\n')

    cat0("  P-Values:\n")
    print(quantile(object[[cft]]$pvalues))

    cat0(statInner("  Cluster Map: ", object[[cft]]$clustermap),
         statInner("  P Map:       ", object[[cft]]$pmap)
    )
  }
}

#' Image a pbj object
#' See image.statMap for additional arguments
#'
#' @export
#' @param x pbj object to create images for
#' @param alpha numeric; threshold to apply to adjusted p-values.
#' @param ... Arguments passed to image.statMap
#' @importFrom graphics image
#' @importFrom graphics par
image.pbj <- function(x, alpha=0.05, ...)
{
  x$mask = if(is.character(x$mask)) readNifti(x$mask) else x$mask
  for(cft in names(x)[ ! names(x) %in% c('stat', 'template', 'mask', 'df') ]){
    stat = x$stat
    # mask stat image with significant voxels
    stat[ abs(x[[cft]]$pmap) < (-log10(alpha) ) ] = 0
    # create a barebones statmap object
    statmap = list(stat=stat, sqrtSigma=NULL, mask=x$mask, template=x$template, formulas=NULL, robust=NULL)
    class(statmap) = "statMap"
    # call image.statMap
    if(x$df==0){
      image(statmap, thresh=qnorm(1-as.numeric(gsub("[^0-9\\.]", "", cft)) ), ...  )
    } else {
      image(statmap, thresh=qchisq(as.numeric(gsub("[^0-9\\.]", "", cft)) * 2, df=df, lower.tail=FALSE ), ...  )
    }
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
  if(!file.exists(outdir)) dir.create(outdir)
  statimg = file.path(outdir, 'stat.nii.gz')
  if(!file.exists(statimg)){
    if(is.character(x$stat)){
      file.copy(x$stat, statimg)
    } else {
      writeNifti(x$stat, statimg)
    }
  }
  if(is.null(x$template)) x$template = x$mask
  sform = do.call(rbind, niftiHeader(x$template)[c('srow_x', 'srow_y', 'srow_z')])
  voxvol = prod(pixdim(x$template))
  for(cft in names(x)[ ! names(x) %in% c('stat', 'template', 'mask', 'df') ]){
    pmapimg = file.path(outdir, paste0('pbj_sei_log10p_', cft, '.nii.gz'))
    clustmapimg = file.path(outdir, paste0('pbj_sei_clust_', cft, '.nii.gz'))
    writeNifti(x[[cft]]$pmap, pmapimg)
    writeNifti(x[[cft]]$clustermap, clustmapimg)

    ### WRITE OUT CLUSTER STATISTICS TABLE ###
    clustmapinds = unique(c(x[[cft]]$clustermap))
    clustmapinds = sort(clustmapinds[ clustmapinds>0])
    clusttab = data.frame('Index'=numeric(0), 'Adjusted p-value'=numeric(0), 'Volume (mm)'=numeric(0), 'Centroid'= character(0), stringsAsFactors = FALSE, check.names = FALSE)
    tabname = file.path(outdir, paste0('sei_table_', cft, '.csv') )
    for(ind in clustmapinds){
      clusttab[ind,c('Index','Adjusted p-value', 'Volume (mm)')] = c(ind, 10^(-x[[cft]]$pmap[ which(x[[cft]]$clustermap==ind) ][1]), sum(x[[cft]]$clustermap==ind)*voxvol)
      clusttab[ind, 'Centroid'] = paste(round(sform %*% c(colMeans(which(x[[cft]]$clustermap==ind, arr.ind=TRUE)), 1 ), 0), collapse=', ')
    }
    write.csv(clusttab, row.names=FALSE, file=tabname)
  }
}

#' Image a CoPE object
#'
#' See image.statMap for additional arguments
#'
#' @export
#' @param x pbj object to create images for
#' @param alpha Displays 1-alpha CoPE maps
#' @param ... Arguments passed to image.statMap
#' @importFrom graphics image
#' @importFrom graphics par
image.CoPE <- function(x, alpha=0.05, ...)
{
  for(rsq in names(x)[ ! names(x) %in% c('stat', 'template') ]){
    statminus = statplus = x$stat
    x[[rsq]]
    # mask stat image with significant voxels
    statminus[ x[[rsq]]$Aminus > 1-alpha ] = 0
    # create a barebones statmap object
    statmap = list(stat=statminus, sqrtSigma=NULL, mask=NULL, template=x$template, formulas=NULL, robust=NULL)
    class(statmap) = "statMap"
    # call image.statMap
    image(statmap, thresh=0.01  )

    # mask stat image with significant voxels
    statplus[ x[[rsq]]$Aplus <= 1-alpha ] = 0
    # create a barebones statmap object
    statmap = list(stat=statplus, sqrtSigma=NULL, mask=NULL, template=x$template, formulas=NULL, robust=NULL)
    class(statmap) = "statMap"
    # call image.statMap
    image(statmap, thresh=0.01  )
  }
}

