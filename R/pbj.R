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
      image(statmap, thresh=qchisq(as.numeric(gsub("[^0-9\\.]", "", cft)) * 2, df=x$df, lower.tail=FALSE ), ...  )
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
#' @importFrom utils write.csv
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
  sform = do.call(rbind, RNifti::niftiHeader(x$template)[c('srow_x', 'srow_y', 'srow_z')])
  voxvol = prod(RNifti::pixdim(x$template))
  for(cft in names(x)[ ! names(x) %in% c('stat', 'template', 'mask', 'df') ]){
    pmapimg = file.path(outdir, paste0('pbj_sei_log10p_', cft, '.nii.gz'))
    clustmapimg = file.path(outdir, paste0('pbj_sei_clust_', cft, '.nii.gz'))
    writeNifti(x[[cft]]$pmap, pmapimg)
    writeNifti(x[[cft]]$clustermap, clustmapimg)

    ### WRITE OUT CLUSTER STATISTICS TABLE ###
    clustmapinds = unique(c(x[[cft]]$clustermap))
    clustmapinds = sort(clustmapinds[ clustmapinds>0])
    clusttab = data.frame('Index'=numeric(0), 'Adjusted p-value'=numeric(0), 'Signed log10(p-value)'=numeric(0), 'Volume (mm)'=numeric(0), 'Centroid'= character(0), stringsAsFactors = FALSE, check.names = FALSE)
    tabname = file.path(outdir, paste0('sei_table_', cft, '.csv') )
    for(ind in clustmapinds){
      clusttab[ind,c('Index','Adjusted p-value', 'Signed log10(p-value)', 'Volume (mm)')] = c(ind, 10^(-abs(x[[cft]]$pmap[ which(x[[cft]]$clustermap==ind) ][1])), x[[cft]]$pmap[ which(x[[cft]]$clustermap==ind) ][1], sum(x[[cft]]$clustermap==ind)*voxvol)
      clusttab[ind, 'Centroid'] = paste(round(sform %*% c(colMeans(which(x[[cft]]$clustermap==ind, arr.ind=TRUE)), 1 ), 0), collapse=', ')
    }
    write.csv(clusttab[order(clusttab[, 'Adjusted p-value']),], row.names=FALSE, file=tabname)
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

#' A vectorized version of papx_edgeworth
#'
#' See image.statMap for additional arguments
#'
#' @param stat Vector of test statistics
#' @param mu3 The third moment of the test statistics
#' @param mu4 The fourth moment of the test statistics
#' @importFrom PDQutils papx_edgeworth
#' @importFrom PDQutils moment2cumulant
vpapx_edgeworth = Vectorize(function (stat, mu3, mu4) PDQutils::papx_edgeworth(stat, raw.cumulants=PDQutils::moment2cumulant(c(0,1, mu3, mu4) ) ))


#' Computes contiguous clusters from a statistical image given a threshold
#'
#' @param stat A statistical Nifti image as an RNifti image object.
#' @param mask A statistical Nifti image mask used in the analysis or a character path to one.
#' @param thr Vector of thresholds to threshold the test statistic image.
#' @param method character string 'extent' or 'mass' indicating whether the cluster extent or cluster mass statistic should be used.
#' @param kernel The kernel type to compute connected components.
#' @param rois If TRUE, return image with voxel values having the indices of the clusters returned if rois=FALSE.
#' @return Returns list of tables of sizes of the connected components above thr.
#' @export
#' @importFrom mmand shapeKernel
#'
cluster = function(stat, mask, thr, method=c('extent', 'mass'), kernel='box', rois=FALSE){
  method = tolower(method[1])
  if(is.character(mask)) mask = readNifti(mask)
  ndims = length(dim(mask))
  tmp = mask
  k = mmand::shapeKernel(3, ndims, type=kernel)
  tmp = lapply(thr, function(th){ tmp[ mask!=0] = (stat[mask!=0]>th); tmp})
  if(rois){
    ccomps = lapply(tmp, function(tm){cc = mmand::components(tm, k); mask[mask!=0]=0; mask[!is.na(cc)] = cc[!is.na(cc)]; mask} )
  } else {
    ccomps = switch(method,
                    'extent'={
                      lapply(tmp, function(tm) table(c(mmand::components(tm, k))) )
                    },
                    'mass'={
                      lapply(tmp, function(tm) c(by(c(stat), c(mmand::components(tm, k)), sum) ))
                    })
    names(ccomps) = paste0(method, '_cft', thr)
  }
  return(ccomps)
}

#' Computes local maxima from an nifti image
#'
#' @param stat A statistical Nifti image as an RNifti image object.
#' @param kernel Type of kernel to use for max/dilation filter
#' @param width Width of kernel (assumes isotropic)
#' @param rois If TRUE, return image with voxel values having the indices of the local maxima returned if rois=FALSE.
#' @return Returns list of tables of sizes of the connected components above thr.
#' @export
#' @importFrom mmand shapeKernel
#'
maxima = function(stat, kernel='box', width=7, rois=FALSE){
  ndims = length(dim(stat))
  dil = dilate(stat, shapeKernel(width, ndims, type=kernel))
    stat[which(stat<dil)] = 0
    imginds = which(stat!=0)
    if(rois){
      stat[imginds] = 1:length(imginds)
      stat
    } else {
      stat[ imginds]
    }
}

#' Computes empirical weighted cdf. Modified from ecdf
#'
#' @param x vector of values
#' @param w vector with length(w)=length(x) of weights
#' @return Returns list of tables of sizes of the connected components above thr.
#' @export
#' @importFrom stats approxfun
#'
wecdf = function (x, w=rep(1, length(x)))
{
  o = order(x)
  x <- x[o]
  w <- w[o]
  n <- length(x)
  sw <- sum(w)
  if (n < 1)
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  # by command sorts x values (again)
  rval <- approxfun(vals, cumsum(c(by(w, x, sum) ))/sw,
                    method = "constant", yleft = 0, yright = 1-1/n, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}


#' Create cluster/maxima summary table
#'
#' @param x pbj object
#' @param method What statistic to provide summary for? must have run that
#' analysis using the pbjInference and mmeStat functions.
#' @return Returns table of unadjusted and FWER adjusted p-values and other summary statistics. Results depend on what statistic
#' function was used for pbjInference.
#' @seealso [mmeStat], [pbjInference]
#' @export
#'
table.pbj = function(x, method=c('CEI', 'maxima', 'CMI')){
  method = method[1]
  ind = grep(method, names(x$obsStat))
  Table = data.frame('Cluster Extent' = c(x$obsStat[[ind]]),
                     'Centroid (vox)' = sapply(1:length(x$obsStat[[ind]]), function(ind) paste(round(colMeans(which(x$ROIs[[2]]==ind, arr.ind = TRUE) )), collapse=', ' )), #RNifti::voxelToWorld(
                     'Unadjusted p-value' = (1-x$margCDF[[ind]](c(x$obsStat[[ind]]))),
                     'FWER p-value' = (1-x$globCDF[[ind]](c(x$obsStat[[ind]]))),
                     check.names=FALSE )
  names(Table) = if(method=='CEI') c('Cluster Extent', 'Centroid (vox)',
'Unadjusted p-value', 'FWER p-value') else if(method=='maxima') c('Chi-square',
'Coord (vox)', 'Unadjusted p-value', 'FWER p-value') else if(method=='CMI')
c('Cluster Mass', 'Centroid (vox)', 'Unadjusted p-value', 'FWER p-value')
  Table = Table[order(Table[,1], decreasing = TRUE),]
  Table
}

#' Compute maxima, CMI, or CEI inference statistics
#'
#' @param stat statistical image
#' @param rois passed to maxima and cluster functions. Return image with ROI
#' indices?
#' @param mask Mask image.
#' @param max Compute local maxima?
#' @param CMI Compute cluster masses?
#' @param CEI Compute cluster extents?
#' @param thr Threshold for CEI or CMI.
#' @return Returns a list with the maxima and CEI for the given image.
#' function was used for pbjInference.
#' @export
#'
mmeStat = function(stat, rois=FALSE, mask, thr, max=FALSE, CMI=FALSE,
CEI=TRUE){
  res = c()
  if(max){
   res = c(maxima=list(maxima(stat, rois=rois)) )
  }
  if(CMI) {
    res = c(res, CMI=cluster(stat, mask=mask, thr=thr, rois=rois, method='mass'))
  }
  if(CEI){
    res = c(res, CEI=cluster(stat, mask=mask, thr=thr, rois=rois,
method='extent'))
  }
  res
}


# color bar function
colorBar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, round(ticks, 2), las=1, cex.axis=cex*0.7, font=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


# image.pbj(obj, slices=NULL, roi=NULL, outputdir=NULL, orientation=c('axial', 'sagittal', 'coronal') ){
#   image(displayParamBoot, template, thresh=(-log10(0.01)), index=slice, cex=cex)
#   if(!is.null(slices)){
#     for (slice in slices){
#       if(!is.null(outputdir)){
#         fname = file.path(outputdir, 'images', paste0('slice', slice, '.png') )
#         dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
#         png(filename = fname, height=4, width=4, units = 'in', res = 300)
#       }
#       cex=1.5
#       # graphical parameters
#       fgcol = 'white'
#       bgcol = 'black'
#       oldpar = par(mgp=c(0.9,.7,0), lwd=1.5, lend=2,
#           cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
#           mar=c(0,2.2,2.2,0), bty='l', oma=c(0,0,0,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)
#       layout(cbind(matrix(1:4, nrow=2, byrow=TRUE) %x% matrix(1, nrow=3, ncol=3) , 5))
#       # display parametric bootstrap statistic
#
#       mtext('Parametric', side=2, cex = 0.8*cex, font = 2)
#       mtext('Bootstrap', side = 3, cex = 0.8*cex, font = 2)
#       # display parametric permutation statistic
#       image(displayParamPerm, template, thresh=(-log10(0.01)), index=slice, cex=cex)
#       mtext('Permutation', side = 3, cex = 0.8*cex, font=2)
#       # display Robust bootstrap statistic
#       image(displayRobustBoot, template, thresh=(-log10(0.01)), index=slice, cex=cex)
#       mtext('Robust', side=2, cex = 0.8*cex, font = 2 )
#       # display robust permutation statistic
#       image(displayRobustPerm, template, thresh=(-log10(0.01)), index=slice, cex=cex)
#
#       # main title
#       #mtext('Probability', side=3, outer = TRUE, cex=1*cex, font=2)
#
#       par(mar=c(10,2,10,0.5))
#       color.bar(pbj:::redyellow(64), min=threshs[1], max=threshs[2], nticks=4)
#       if(!is.null(outputdir)) dev.off()
#     }
#   }
# }


