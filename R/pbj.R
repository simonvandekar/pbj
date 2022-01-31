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
#' @param cft A a vector of cluster forming thresholds for the test statistic image. Will compute cluster sizes for each threshold.
#' @param method character string 'extent' or 'mass' indicating whether the cluster extent or cluster mass statistic should be used.
#' @param kernel The kernel type to compute connected components.
#' @param rois If TRUE, return image with voxel values having the indices of the clusters returned if rois=FALSE.
#' @return Returns list of tables of sizes of the connected components above thr.
#' @export
#' @importFrom mmand shapeKernel
#'
cluster = function(stat, mask, cft, method=c('extent', 'mass'), kernel='box', rois=FALSE){
  method = tolower(method[1])
  if(is.character(mask)) mask = readNifti(mask)
  ndims = length(dim(mask))
  tmp = mask
  k = mmand::shapeKernel(3, ndims, type=kernel)
  tmp = lapply(cft, function(th){ tmp[ mask!=0] = (stat[mask!=0]>th); tmp})
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
    # modifies ccomps attribute in cluster function
    ccomps = lapply(1:length(cft), function(ind){ attributes(ccomps[[ind]]) <- list('cft'=cft[ind]); ccomps[[ind]]})
  }
  return(ccomps)
}

#' Computes local maxima from an nifti image
#'
#' @param stat A statistical Nifti image as an RNifti image object.
#' @param kernel Type of kernel to use for max/dilation filter
#' @param width Width of kernel (assumes isotropic)
#' @param rois If TRUE, return image with voxel values having the indices of the local maxima returned if rois=FALSE.
#' @return Returns vector of local maxima in the image.
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
#' @param cft cluster forming threshold to display. If NULL, just display the first.
#' @return Returns table of unadjusted and FWER adjusted p-values and other summary statistics. Results depend on what statistic
#' function was used for pbjInference.
#' @seealso [mmeStat], [cluster], [maxima], [pbjInference]
#' @export
#'
table.statMap = function(x, method=c('CEI', 'maxima', 'CMI'), cft=NULL){
  x = x$pbj
  if(is.null(x)){
    stop("Must run pbjInference to create cluster table.")
  }
  method = method[1]
  ind = grep(method, names(x$obsStat))
  # get indices corresponding to this method
  cfts = sapply(x$obsStat[ind], attr, which='cft')
  # maxima don't have cft attribute
  if(!is.null(cft) & method!='maxima'){
    ind = ind[which(cfts==cft)]
    if(length(ind)==0){
      stop('Specified cft is ', cft, '. Existing cfts are ', paste(cfts, collapse=', '))
    }
  }
  ind = ind[1]
  Table = data.frame('Cluster Extent' = c(x$obsStat[[ind]]),
                     'Centroid (vox)' = sapply(1:length(x$obsStat[[ind]]), function(ind2) paste(round(colMeans(which(x$ROIs[[ind]]==ind2, arr.ind = TRUE) )), collapse=', ' )), #RNifti::voxelToWorld(
                     'Unadjusted p-value' = (1-x$margCDF[[ind]](c(x$obsStat[[ind]]))),
                     'FWER p-value' = (1-x$globCDF[[ind]](c(x$obsStat[[ind]]))),
                     check.names=FALSE )
  names(Table) = if(method=='CEI') c('Cluster Extent', 'Centroid (vox)',
                                     'Unadjusted p-value', 'FWER p-value') else if(method=='maxima') c('Chi-square',
                                                                                                       'Coord (vox)', 'Unadjusted p-value', 'FWER p-value') else if(method=='CMI')
                                                                                                         c('Cluster Mass', 'Centroid (vox)', 'Unadjusted p-value', 'FWER p-value')
  Table = Table[order(Table[,1], decreasing = TRUE),]
  Table[,'cluster ID'] = 1:nrow(Table)
  Table = Table[, c(ncol(Table), 1:(ncol(Table)-1) )]
  #attributes(Table) <- list('cft'=cfts[ind], 'ind'=ind)
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
#' @param cft A single threshold for CEI or CMI.
#' @return Returns a list with the maxima and CEI for the given image.
#' This function is used  as the default `statistic` argument in [pbjInference()].
#' @export
#'
mmeStat = function(stat, rois=FALSE, mask, cft, max=FALSE, CMI=FALSE,
                   CEI=TRUE){
  res = c()
  if(max){
    res = c(maxima=list(maxima(stat, rois=rois)) )
  }
  if(CMI) {
    res = c(res, CMI=cluster(stat, mask=mask, cft=cft, rois=rois, method='mass'))
  }
  if(CEI){
    res = c(res, CEI=cluster(stat, mask=mask, cft=cft, rois=rois,
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



#' View pbj results
#'
#' Uses a pbj object and statMap object to visualize CEI, CMI, or maxima results
#'
#' @param pbjObj pbj object to visualize
#' @param object the statMap that created the pbj object.
#' @param cft The cluster forming threshold or threshold for visualizing results (in the case of maxima). On the chi-square scale.
#' @param pCFT For convenience, the user can specify the CFT in terms of a p-value.
#' @param clusterMask logical indicating whether to mask the results with significant clusters.
#' @param alpha Adjusted p-value threshold to display clusters.
#' @param ... Arguments passed to image.niftiImage.
#' @importFrom utils write.csv
#' @export
image.statMap = function(object, method=c('CEI', 'maxima', 'CMI'), cft=NULL, pCFT=NULL, roi=NULL, slice=NULL, alpha=NULL, clusterMask=TRUE, clusterID=TRUE, outputdir=NULL, plane=c('axial', 'sagittal', 'coronal'), ... ){

  if(!is.null(pCFT)) cft = qchisq(pCFT, df=object$sqrtSigma$df, lower.tail=FALSE)
  # use mask if user didn't provide a template
  if(is.null(object$template)) object$template=object$mask
  # read template if stored as a character
  x = if(is.character(object$template)) readNifti(object$template) else object$template
  method = method[1]
  plane=plane[1]

  pbjObj = object$pbj
  # if user hasn't run pbj yet visualize using image.niftiImage
  if(is.null(pbjObj)){
    if(is.null(cft)) cft=0
    image.niftiImage(stat.statMap(object), bgimg = object$template, thresh = cft, index = slice, plane=plane)
  } else {
    if(is.character(object$mask)){
      object$mask = readNifti(object$mask)
    }
    ind = grep(method, names(pbjObj$obsStat))
    # get indices corresponding to this method
    cfts = sapply(pbjObj$obsStat[ind], attr, which='cft')
    # maxima don't have cft attribute
    if(!is.null(cft) & method!='maxima'){
      ind = ind[which(cfts==cft)]
      if(length(ind)==0){
        stop('Specified cft is ', cft, '. Existing cfts are ', paste(cfts, collapse=', '))
      }
    } else {
      cft = cfts[1]
    }
    ind = ind[1]
    if(is.null(cft)) stop("Must specify cft for visualization with method='maxima'")

    st = table.pbj(pbjObj, method=method, cft=cft)
    imgdims = dim(stat.statMap(object))
    planenum = switch(plane, "axial"=3, 'sagittal'=1, 'coronal'=2)
    otherplanes = which(!1:3 %in% planenum)
    if(!is.null(alpha)){
      pbjObj$ROIs[[ind]][ pbjObj$ROIs[[ind]] %in% which(st$`FWER p-value`>alpha) ] = 0
      object$stat[pbjObj$ROIs[[ind]][ object$mask>0 ] ==0] = 0
      st = st[which(st$`FWER p-value`<=alpha),]
    } else if(clusterMask){
      # zero values outside of clusterMask
      object$stat[pbjObj$ROIs[[ind]][ object$mask>0 ]==0] = 0
    }

    # image.statMap crops the image, so we need to do that here
    xoffset = which(apply(x!=0, 1, any))
    xoffset = c(xoffset[1], length(xoffset))+ c(-1, 0)
    yoffset = which(apply(x!=0, 2, any))
    yoffset = c(yoffset[1], length(yoffset))+ c(-1, 0)
    zoffset = which(apply(x!=0, 3, any))
    zoffset = c(zoffset[1], length(zoffset))+ c(-1, 0)
    offsets = list(xoffset, yoffset, zoffset)
    rm(zoffset, yoffset, xoffset)
    # if roi is non-null draw single slice with centroid of each selected ROI
    if(!is.null(roi)){
      for(ro in roi){
        # get coordinates from table
        coords = as.numeric(strsplit(st[ro,3], split=', ')[[1]])

        # slice to display
        planecoord = coords[planenum]
        # coordinates to display the cluster index number
        othercoords = coords[-planenum]
        if(clusterID){
          otherfunc=function(){text(othercoords[1]-offsets[[otherplanes[1]]][1]+1, othercoords[2]-offsets[[otherplanes[2]]][1]+1, labels=ro, col='white', font=2)}
          # trying to debug coordinates
          #otherfunc=function(){points(othercoords[1]-offsets[[otherplanes[1]]][1], othercoords[2], col='white')}
          #otherfunc = function(){text(50, seq(1, offsets[[otherplanes[2]]][2], 10), labels=seq(1, offsets[[otherplanes[2]]][2], 10), col='white')}
          #otherfunc = function(){text(seq(1, offsets[[otherplanes[1]]][2], 10), 50, labels=seq(1, offsets[[otherplanes[1]]][2], 10), col='white')}
          image(stat.statMap(object), bgimg=object$template, plane=plane, index=coords[planenum], thresh=cft, other=otherfunc, ...)
        } else {
          image(stat.statMap(object), bgimg=object$template, plane=plane, index=coords[planenum], thresh=cft, ...)
        }
      }
    } else if(!is.null(slice)){
      for (slic in slice){
        if(clusterID){
          coords = do.call(rbind, lapply(strsplit(st[,3], split=', '), as.numeric))
          coordInds = which(coords[,planenum]==slic)
          if(length(coordInds)==0){
            image(stat.statMap(object), bgimg=object$template, plane=plane, index=slic, thresh=cft, ...)
          } else {
            coordLabels = st$`cluster ID`[coordInds]
            othercoords = coords[coordInds,-planenum, drop=FALSE]
            # if(plane %in% c('axial', 'coronal'){
            #   #p1 =
            #   #p2 =
            # } else {
            #
            # }
            otherfunc=function(){text(othercoords[,1]-offsets[[otherplanes[1]]][1], othercoords[,2]-offsets[[otherplanes[2]]][1], labels=coordLabels, col='white', font=2)}
            image(stat.statMap(object), bgimg=object$template, plane=plane, index=slic-offsets[[planenum]][1], thresh=cft,
                  other=otherfunc, ...)
          }
        } else {
          image(stat.statMap(object), bgimg=object$template, plane=plane, index=slic-offsets[[planenum]][1], thresh=cft, ...)
        }
      }

      # if roi and slice are null, do lightbox view
    } else {
      image(stat.statMap(object), bgimg=object$template, plane=plane, thresh=cft, ...)

    }
  }
}


