#' @importFrom utils str
#' @export
summary.statMap <- function(object, ...)
{
  cat("Full formula: ", paste0(as.character(object$formulas[[1]])) )
  cat("\nReduced formula: ", as.character(object$formulas[[2]]), collapse='', '\n')
  cat(sum(if(is.character(object$mask)) readNifti(object$mask) else object$mask), ' voxels in mask')
  cat("\nStatMap quantiles: [", paste(round(quantile(object$stat, probs = c(0, 0.01, 0.05, 0.95, 0.99, 1)), 2), collapse=', '), "]" )
  cat("\nsqrtSigma: \n")
  cat("  [n = ", object$sqrtSigma$n, '; df = ', object$sqrtSigma$df, '; rdf = ', object$sqrtSigma$rdf, ']\n')
  cat("id variable is:\n")
  str(object$sqrtSigma$id)
  cat("lmPBJ inference settings:\n  robust = ", object$sqrtSigma$robust, '; HC3 = ', object$sqrtSigma$HC3, '; transform = ',  object$sqrtSigma$transform)
  if(is.null(object$pbj)){
    cat('\npbjInference not run yet.')
  } else {
    cat('\npbjInference results:')
    if(any(grepl('CEI', names(object$pbj$obsStat))) ){
      cat('\n  CEI run with cft = ', paste(round(sapply(object$pbj$obsStat[grep('CEI', names(object$pbj$obsStat))], attr, 'cft'), 2), collapse=', ') )
    }
    if(any(grepl('CMI', names(object$pbj$obsStat)) )){
      cat('\n  CMI run with cft = ', paste(round(sapply(object$pbj$obsStat[grep('CMI', names(object$pbj$obsStat))], attr, 'cft'), 2), collapse=', ') )
    }
    if(any(grepl('maxima', names(object$pbj$obsStat)))){
      cat('\n  Inference run on local maxima.')
    }
  }
}

#' @export
#' @method print statMap
print.statMap <- function(x, ...)
{
  summary(x, ...)
}

redyellow = colorRampPalette(c('red', 'yellow'), space='Lab')
bluecyan = colorRampPalette(c('blue', 'cyan'), space='Lab')


#' Write the statMap objects out
#'
#' Given a statMap object and a directory write the objects as stat.nii.gz, coef.nii.gz and sqrtSigma.nii.gz
#' @param x the statMap object to write out
#' @param outdir the directory to write into
#' @param images Logical indicating whether or not to write out the nifti images.
#' @param sqrtSigma Logical indicating whether or not to write out the objects needed for bootstrapping.
#' @return a list of what was written
#' @export
write.statMap <- function(x, outdir, images=TRUE, sqrtSigma=TRUE){
  statimg  = file.path(outdir, 'stat.nii.gz')
  coefimg   = file.path(outdir, 'coef.nii.gz')
  res   = file.path(outdir, 'sqrtSigma.rds')
  if(is.character(x$stat)){
    if(images){
      file.copy(x$stat, statimg)
      file.copy(x$sqrtSigma, res)
      file.copy(x$coef, coefimg)
    }
  } else {
    if(images){
      message('Writing stat and coef images.')
      dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
      writeNifti(stat.statMap(x), statimg)
      writeNifti(coef.statMap(x), coefimg)
    }

    if(sqrtSigma){
      message('Writing sqrtSigma object.')
      saveRDS(x$sqrtSigma, file = res)
    }
  }


  inference=list()
  pbj = x$pbj
  if(!is.null(pbj)){
    #if(is.null(x$template)) x$template = x$mask
    #sform = do.call(rbind, RNifti::niftiHeader(x$template)[c('srow_x', 'srow_y', 'srow_z')])
    #voxvol = prod(RNifti::pixdim(x$template))
    message('Writing inference images.')
    for(infType in names(pbj$ROIs)){
      if(grepl('CEI|CMI', infType)) cft = attr(pbj$obsStat[[infType]], 'cft')
      # removes number from name if infType
      method = gsub("[0-9]", '', infType)
      tab = table.statMap(x, method=method, cft=cft)
      pmapimg = file.path(outdir, paste0('log10p_', infType, '.nii.gz'))
      clustmapimg = file.path(outdir, paste0('clustIDs_', infType, '.nii.gz'))
      pmap = pbj$ROIs[[infType]]
      # sets clusters to their p-value. This is wrong currently
      pmap[ which(pmap!=0) ] = -log10(tab$`FWER p-value`[ pmap[which(pmap!=0)] ])
      writeNifti(pmap, pmapimg)
      writeNifti(pbj$ROIs[[infType]], clustmapimg)
      inference[[infType]] = c(pmapimg, clustmapimg)
    }
  }
  return(c(list(stat=statimg, coef=coefimg, sqrtSigma=res), inference) )
}

#' Gets a 4D niftiImage of the coefficient image from a statMap object
#'
#' Returns a statistical niftiImage object from a statMap object
#' @param x the statMap object to extract a coefficient niftiImage from
#' @return a niftiImage object of the coefficient image
#' @export
stat.statMap = function(x){
  if(is.character(x$stat)){
    stat = readNifti(x$stat)
  } else {
    # output 4D coefficient image
    stat = x$mask
    stat[ stat!=0 ] = x$stat
  }
  return(stat)
}

#' Gets a 4D niftiImage of the statistical image from a statMap object
#'
#' Returns a 4D coefficient niftiImage object from a statMap object
#' @param object the statMap object to extract a coefficient niftiImage from
#' @param ... additional arguments (ignored)
#' @return a niftiImage object of the coefficient image
#' @export
coef.statMap = function(object, ...){
  # output 4D coefficient image
  coef = simplify2array(lapply(1:nrow(object$coef), function(coefv){ object$mask[object$mask!=0] = object$coef[coefv,]; object$mask}))
  coef = updateNifti(coef, template=object$mask)
  return(coef)
}

#' Gets a niftiImage of the variance image from a statMap object
#'
#' Will return a niftiImage of the variance image
#' @param x the statMap to extract the variance image from
#' @return a niftiImage object of the variance image
#' @export
var.statMap = function(x){
  if(x$df>1)
    stop('Only supported for df<1.')
  # get niftiImage from mask
  varimg = x$mask
  varimg[ varimg!=0] = rowSums(x$sqrtSigma$res^2)/x$sqrtSigma$rdf
  return(varimg)
}

#' Plots the tested variables in a statMap object
#'
#' Returns a statistical niftiImage object from a statMap object
#' @param x the statMap object with pbjInference run
#' @param emForm a formula specifying how to plot the data
#' @param method A character specifying which inference method to plot results for. One of maxima, CEI, or CMI
#' @param cft A numeric specifying CFT to plot results for, defaults to the first one
#' @param roiInds A numeric/integer vector specifying which ROIs to plot results for.
#' @param ... arguments passed to plot
#' @return a niftiImage object of the coefficient image
#' @importFrom graphics polygon points
#' @importFrom stats terms
#' @importFrom emmeans emmeans ref_grid
#' @importFrom scales alpha
#' @export
plot.statMap = function(x, emForm=NULL, method='CEI', cft=NULL, roiInds=NULL, ...){

  pbjObj = x$pbj
  if(is.null(pbjObj)) stop('Must run pbjInference to plot results.')
  st = table.statMap(x, method, cft)
  ind = inferenceIndex(x$pbj$obsStat, method=method, cft=cft)
  rois = pbjObj$ROIs[[ind]]
  obsstat = pbjObj$obsStat[[ind]]


  # load data
  imgs = RNifti::readNifti(x$images)
  # fit model on average data
  full = as.formula(x$formulas$full)
  red = as.formula(x$formulas$reduced)
  fullT = terms(full)
  redT = terms(red)
  # formula elements to condition on
  term = attr(fullT, 'term.labels')[!attr(fullT, 'term.labels') %in% attr(redT, 'term.labels') ]
  condVars = sapply(all.vars(full), function(x) grepl(x, term) )
  condVars = names(condVars)[condVars]
  robust = x$sqrtSigma$robust
  data = get_all_vars(full, data=data)
  for(roiInd in roiInds){
    data$y = sapply(imgs, function(img) mean(img[ rois==roiInd]) )

    fullModel = lm(update.formula(full, y ~ .), data=data)
    #redModel = lm(update.formula(red, y ~ .), data=data)

    # emmeans argument using condVars
    # known warning about as.numeric on characters
    suppressWarnings(ats <- apply(data[,condVars], 2, function(x) sort(unique(ifelse(is.na(as.numeric(x)), x, as.numeric(x) ) ) ) ))
    emg = ref_grid(fullModel, at=ats)
    plotdf = as.data.frame(emmeans(emg, ~ age | sex))


    cex=1.5
    # graphical parameters
    fgcol = 'white'
    bgcol = 'black'
    par(mgp=c(1.7,.7,0), lwd=1.5, lend=2,
        cex.lab=cex, cex.axis=0.8*cex, cex.main=1*cex,
        mar=c(3,3,2.2,0), bty='l', oma=c(0,0,0,0), bg=bgcol, fg=fgcol, col.axis=fgcol, col.lab=fgcol, col.main = fgcol, col.sub=fgcol)

    # plot predictions with raw data
    cols=c('#fb9a99', '#a6cee3', '#e31a1c', '#1f78b4')
    plot(data[,'age'], data[,'y'], col=cols[2+as.numeric(factor(data[,'sex']))], pch=16, ylab='Gray matter volume (AU)', xlab='Age (Years)', ...)
    by(plotdf, plotdf$sex, function(df){
      polygon(c(df$age, rev(df$age)), c(df$lower.CL, rev(df$upper.CL)), col=alpha(cols[as.numeric(df$sex)], alpha=0.8), border=NA)
      points(df$age, df$emmean, type='l', col=cols[as.numeric(df$sex)+2])
    } )
  }

}



#' Plots the tested variables in a statMap object
#'
#' Returns a statistical niftiImage object from a statMap object
#' @param statMap the statMap object with pbjInference run
#' @param method A character specifying which inference method to plot results for. One of maxima, CEI, or CMI
#' @param cft A numeric specifying CFT to plot results for, defaults to the first one
#' @param roiInds A numeric/integer vector specifying which ROIs to plot results for.
#' @param data Data frame with covariates. Optional if available in statMap object.
#' @return a niftiImage object of the coefficient image
#' @export
roiMeans = function(statMap, method='CEI', cft=NULL, roiInds=NULL, data=NULL){
  pbjObj = statMap$pbj
  full = statMap$formulas$full
  if(is.null(pbjObj)) stop('Must run pbjInference to plot results.')
  st = table.statMap(statMap, method, cft)
  ind = inferenceIndex(statMap$pbj$obsStat, method=method, cft=cft)
  rois = pbjObj$ROIs[[ind]]
  obsstat = pbjObj$obsStat[[ind]]

  # load data
  imgs = RNifti::readNifti(statMap$images)
  if(is.null(data)){
    data = get_all_vars(full, statMap$data)
  } else {
    data = get_all_vars(full, data=data)
  }
  for(roi in roiInds){
    data[,paste0('roi', roi)] = sapply(imgs, function(img) mean(img[ rois==roi]) )
    return(data)
  }
}

#' Draws a color bar in a figure
#'
#' Can be used in combination with layout and image.statMap to create figures.
#'
#' Returns a statistical niftiImage object from a statMap object
#' @param lut the range of colors to display
#' @param min Where to put the minimum tick mark.
#' @param max Where to put the maximum tick mark.
#' @param nticks Number of ticks on the color bar.
#' @param ticks Locations for ticks.
#' @param ylab label for the color bar.
#' @param las argument passed to par.
#' @param title The title for the color bar
#' @return a niftiImage object of the coefficient image
#' @importFrom graphics axis rect
#' @export
colorBar <- function(lut, min, max=-min, nticks=4, ticks=seq(min, max, len=nticks), title='', ylab='', las=1) {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, round(ticks, 2), las=las, font=2)
  title(ylab=ylab, font=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

#' View statMap and pbj results in a statMap object
#'
#' Uses a statMap with inferences results to visualize CEI, CMI, or maxima
#'
#' @param x the statMap that created the pbj object.
#' @param method Which inference method to visualize.
#' @param cft The cluster forming threshold or threshold for visualizing results (in the case of maxima). On the chi-square scale.
#' @param pCFT For convenience, the user can specify the CFT in terms of a p-value.
#' @param roi Which ROI index to visualize.
#' @param slice Which slice number to visualize.
#' @param clusterID logical indicating whether to include the clusterID in the figure.
#' @param plane Which anatomical plane to visualize.
#' @param clusterMask logical indicating whether to mask the results with significant clusters.
#' @param alpha Adjusted p-value threshold to display clusters.
#' @param oma argument passed to par
#' @param mar argument passed to par
#' @param bg argument passed to par
#' @param ... Arguments passed to image.niftiImage.
#' @importFrom utils write.csv
#' @importFrom graphics text
#' @export
image.statMap = function(x, method=c('CEI', 'maxima', 'CMI'), cft=NULL, pCFT=NULL, roi=NULL, slice=NULL, alpha=NULL, clusterMask=TRUE, clusterID=TRUE, plane=c('axial', 'sagittal', 'coronal'), oma = rep(0, 4), mar = rep(0, 4), bg = "black", ... ){
  if(!is.null(pCFT)) cft = qchisq(pCFT, df=x$sqrtSigma$df, lower.tail=FALSE)
  # set graphical parameters
  par(oma = oma,
      mar = mar, bg = bg)
  oldpar <- par(no.readonly = TRUE)
  # use mask if user didn't provide a template
  if(is.null(x$template)) x$template=x$mask
  # read template if stored as a character
  x$template = if(is.character(x$template)) readNifti(x$template) else x$template
  method = method[1]
  plane=plane[1]

  pbjObj = x$pbj
  # if user hasn't run pbj yet visualize using image.niftiImage
  if(is.null(pbjObj)){
    if(is.null(cft)) cft=0
    image.niftiImage(stat.statMap(x), bgimg = x$template, thresh = cft, index = slice, plane=plane)
  } else {
    if(is.character(x$mask)){
      x$mask = readNifti(x$mask)
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

    st = table.statMap(x, method=method, cft=cft)
    imgdims = dim(stat.statMap(x))
    planenum = switch(plane, "axial"=3, 'sagittal'=1, 'coronal'=2)
    otherplanes = which(!1:3 %in% planenum)
    if(!is.null(alpha)){
      pbjObj$ROIs[[ind]][ pbjObj$ROIs[[ind]] %in% which(st$`FWER p-value`>alpha) ] = 0
      x$stat[pbjObj$ROIs[[ind]][ x$mask>0 ] ==0] = 0
      st = st[which(st$`FWER p-value`<=alpha),]
    } else if(clusterMask){
      # zero values outside of clusterMask
      x$stat[pbjObj$ROIs[[ind]][ x$mask>0 ]==0] = 0
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
          image(stat.statMap(x), bgimg=x$template, plane=plane, index=coords[planenum]-offsets[[planenum]][1], thresh=cft, other=otherfunc, ...)
        } else {
          image(stat.statMap(x), bgimg=x$template, plane=plane, index=coords[planenum]-offsets[[planenum]][1], thresh=cft, ...)
        }
      }
    } else if(!is.null(slice)){
      for (slic in slice){
        if(clusterID){
          coords = do.call(rbind, lapply(strsplit(st[,3], split=', '), as.numeric))
          coordInds = which(coords[,planenum]==slic)
          if(length(coordInds)==0){
            image(stat.statMap(x), bgimg=x$template, plane=plane, index=slic-offsets[[planenum]][1], thresh=cft, ...)
          } else {
            coordLabels = st$`cluster ID`[coordInds]
            othercoords = coords[coordInds,-planenum, drop=FALSE]
            otherfunc=function(){text(othercoords[,1]-offsets[[otherplanes[1]]][1], othercoords[,2]-offsets[[otherplanes[2]]][1], labels=coordLabels, col='white', font=2)}
            image(stat.statMap(x), bgimg=x$template, plane=plane, index=slic-offsets[[planenum]][1], thresh=cft,
                  other=otherfunc, ...)
          }
        } else {
          image(stat.statMap(x), bgimg=x$template, plane=plane, index=slic-offsets[[planenum]][1], thresh=cft, ...)
        }
      }

      # if roi and slice are null, do lightbox view
    } else {
      image(stat.statMap(x), bgimg=x$template, plane=plane, thresh=cft, ...)

    }
  }
}

#' Compute pTFCE image
#'
#' Uses a statMap object to perform pTFCE inference
#'
#' @param statMap A statmap object.
#' @param ... Arguments passed to ptfce.
#' @importFrom pTFCE ptfce
#' @importFrom stats pchisq
#' @importFrom utils capture.output
#' @importFrom oro.nifti as.nifti
#' @export
ptfce.statmap = function(statMap, ...){
  if(is.character(statMap$mask))  statMap$mask = readNifti(statMap$mask)
  statimg = stat.statMap(statMap)
  # convert Chisq to Z
  statimg[statMap$mask>0] = qnorm(pchisq(statimg[statMap$mask>0], df=statMap$sqrtSigma$df, log.p = TRUE), log.p=TRUE)
  invisible(capture.output(test <- ptfce(as.nifti(array(statimg, dim = dim(statimg))), as.nifti(array(statMap$mask, dim=dim(statMap$mask))) )))
  return(test)
}
