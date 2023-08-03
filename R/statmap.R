#' object A statMap object
#' @importFrom utils str
#' @param object A statMap object
#' @param ... additional arguments
#' @export
summary.statMap <- function(object, ...)
{
  cat("Full formula: ", paste0(as.character(object$formulas[[1]])) )
  cat("\nReduced formula: ", as.character(object$formulas[[2]]), collapse='', '\n')
  cat(sum(if(is.character(object$mask)) readNifti(object$mask) else object$mask), ' voxels in mask')
  if(is.character(object$stat)){
    cat("\nStatMap location: \n", object$stat)
  } else {
    cat("\nStatMap quantiles (0, 0.01, 0.05, 0.95, 0.99, 1):\n [", paste(round(quantile(object$stat, probs = c(0, 0.01, 0.05, 0.95, 0.99, 1)), 2), collapse=', '), "]" )
  }
  if(is.character(object$sqrtSigma)){
    cat("\nsqrtSigma location: \n", object$sqrtSigma)
  } else {
    cat("\nsqrtSigma: \n")
    cat("  [n = ", object$sqrtSigma$n, '; df = ', object$sqrtSigma$df, '; rdf = ', object$sqrtSigma$rdf, ']\n')
    cat("id variable is:\n")
    str(object$sqrtSigma$id)
    cat("lmPBJ inference settings:\n  robust = ", object$sqrtSigma$robust, '; HC3 = ', object$sqrtSigma$HC3, '; transform = ',  object$sqrtSigma$transform)
  }
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
#' @param statMethod Passed to stat.statMap to specify how to write the statistical image. p= -log10(p) * sign(coef), S=RESI, chisq=chi-squared statistic
#' @return a list of what was written
#' @export
write.statMap <- function(x, outdir, images=TRUE, sqrtSigma=TRUE, statMethod=c('p', 'S', 'chisq')){
  statMethod = tolower(statMethod[1])
  statimg  = file.path(outdir, paste0('stat_', statMethod, '.nii.gz') )
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
      writeNifti(stat.statMap(x, method = statMethod), statimg)
      writeNifti(coef.statMap(x), coefimg)
    }

    if(sqrtSigma){
      message('Writing sqrtSigma object.')
      saveRDS(x$sqrtSigma, file = res)
    }
  }


  inference=list()
  pbj = x$pbj
  ## SNV: NEED TO EVALUATE THIS OUTPUT
  if(!is.null(pbj)){
    #if(is.null(x$template)) x$template = x$mask
    #sform = do.call(rbind, RNifti::niftiHeader(x$template)[c('srow_x', 'srow_y', 'srow_z')])
    #voxvol = prod(RNifti::pixdim(x$template))
    message('Writing inference images.')
    cft = NULL
    for(infType in names(pbj$ROIs)){
      if(grepl('CEI|CMI', infType)) cft = attr(pbj$obsStat[[infType]], 'cft') else cft=NULL
      # removes number from name if infType
      method = tolower(gsub("[0-9]", '', infType))
      tab = table.statMap(x, method=method, cft_chisq=cft)
      pmapimg = file.path(outdir, paste0('log10p_', gsub(' ', '_', infType), '.nii.gz'))
      clustmapimg = file.path(outdir, paste0('clustIDs_', gsub(' ', '_', infType), '.nii.gz'))
      pmap = pbj$ROIs[[infType]]
      # sets clusters to their p-value.
      for(ind in tab[,1]){
        pmap[ which(pmap==ind) ] = -log10(tab$`FWER p-value`[ tab[,1]==ind ])
      }
      writeNifti(pmap, pmapimg)
      writeNifti(pbj$ROIs[[infType]], clustmapimg)
      inference[[infType]] = c(pmapimg, clustmapimg)
    }
  }
  return(c(list(stat=statimg, coef=coefimg, sqrtSigma=res), inference) )
}

#' Returns a statistical niftiImage object from a statMap object
#'
#'
#' @param x the statMap object to extract a coefficient niftiImage from
#' @param method the scale of visualization, chi-squared statistic, effect size (S=RESI), p-value. If df is equal to 1, the maps are scaled by the sign of the coefficient.
#' @return a niftiImage object of the chi-square statistical image.
#' @importFrom RESI chisq2S
#' @export
#' @examples
#' # loading example data
#' library(pain21)
#' pain = pain21()
#' pdata = pain$data
#'
#' # fitting regression of images onto study sample size, weights proportional to study sample size
#' pbjModel2 = lmPBJ(images=pdata$images, form=~n, formred=~1, W = pdata$n, mask=pain$mask, data=pdata)
#' stat.statMap(pbjModel2)
#' stat.statMap(pbjModel2, method='chisq')
#' stat.statMap(pbjModel2, method='S')
#'
stat.statMap = function(x, method=c('p', 'S', 'chisq')){
  method = tolower(method[1])
  if(is.character(x$stat)){
    stat = readNifti(x$stat)
    res = stat[ x$mask!=0 ]
  } else {
    # output 4D coefficient image
    stat = x$mask
    res = x$stat
    if(method=='s'){
      res = chisq2S(res, x$sqrtSigma$df, x$sqrtSigma$n)
      if(x$sqrtSigma$df==1){
        res = res * sign(x$coef)
      }
    }
    if(method == 'p'){
      res = -pchisq(res, df = x$sqrtSigma$df, lower.tail=FALSE, log.p=TRUE)/log(10)
      if(x$sqrtSigma$df==1){
        res = res * sign(x$coef)
      }
    }
    stat[ stat!=0] = res
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

#' Creates dataframes, emmeans, and models of the selected ROI data
#'
#' Returns a statistical niftiImage object from a statMap object
#' @param x the statMap object
#' @param emForm a formula specifying how to plot the data
#' @param method A character specifying which inference method to plot results for. One of maxima, CEI, or CMI
#' @param cft_s cluster forming threshold on effect size scale.
#' @param cft_p cluster forming threshold on p-value scale.
#' @param cft_chisq cluster forming threshold on chi-square statistic scale.
#' @param roiInds A numeric/integer vector specifying which ROIs to plot results for.
#' @param ... arguments passed to ref_grid
#' @return a list with a dataframe containing the original data with a column added for
#' mean ROI values for each subject, an emmeans objects based on the emForm argument,
#' and model objects of the full model fit to each ROI.
#' @details Returns ROI summary results if pbjInference has been executed on the statMap object.
#' If not, it attempts to create a cluster summary table based on the method and cft_* arguments.
#' emmeans objects can be used for plotting and are returned if emForm is specified.
#' @importFrom emmeans emmeans ref_grid
#' @importFrom RNifti readNifti
#' @importFrom sandwich vcovHC
#' @export
plotData.statMap = function(x, emForm=NULL, method='CEI', cft_s=NULL, cft_p=NULL, cft_chisq=NULL, roiInds=NULL, ...){
  method = tolower(method[1])
  # CFT passed as p value or effect size converted to chi-squared threshold
  cft = cft_chisq
  if(!is.null(cft_p)){
    cft = qchisq(cft_p, df=x$sqrtSigma$df, lower.tail=FALSE)
  }
  if(!is.null(cft_s)){
    cft = cft_s^2 * x$sqrtSigma$n + x$sqrtSigma$df
  }

  data = x$data
  robust = x$sqrtSigma$robust
  data$W = x$sqrtSigma$W
  if(is.null(x$pbj)){
    x$pbj = list()
    mmeStatArgs = as.list(grepl(method, c('maxima', 'cei', 'cmi')))
    names(mmeStatArgs) = c('maxima', 'CEI', 'CMI')
    x$pbj$ROIs = do.call(mmeStat, c(list(stat=stat.statMap(x), rois=TRUE, mask=x$mask, cft=cft), mmeStatArgs) )
    x$pbj$obsStat = do.call(mmeStat, c(list(stat=stat.statMap(x), mask=x$mask, cft=cft), mmeStatArgs) )
  }
  ind = inferenceIndex(x$pbj$obsStat, method=method, cft=cft)
  rois = x$pbj$ROIs[[ind]]
  st = table.statMap(x, method, cft_chisq=cft)

  # load data
  imgs = readNifti(x$images)
  # fit model on average data
  full = as.formula(x$formulas$full)
  red = as.formula(x$formulas$reduced)
  if(is.null(roiInds)){
    roiInds = st[,1]
  }

  data[, paste0('roi', roiInds)] = do.call(cbind,
                                           lapply(roiInds, function(roiInd){
                                             sapply(imgs, function(img) mean(img[ rois==roiInd]) )
                                           }))
  # W is found from data$W
  plotDF = lapply(roiInds, function(roiInd){
    if(is.null(x$sqrtSigma$id)){
      fullModel = lm(update.formula(full, as.formula(paste0('roi', roiInd, '~ .'))), weights=W, data=data)
      #redModel = lm(update.formula(red, y ~ .), data=data)
    }

    if(!is.null(emForm)){
      if(robust){
        # assumes HC3 vcov
        emGrid = ref_grid(fullModel, vcov.=vcovHC, ...)
      } else {
        emGrid = ref_grid(fullModel, ...)
      }
      plotdf = emmeans(emGrid, as.formula(emForm))
    } else {
      plotdf = NULL
    }
    list(plotDF=plotdf, model=fullModel)
  })
  mods = lapply(plotDF, function(x) x$model)
  plotDF = lapply(plotDF, function(x) x$plotDF)
  names(mods) = names(plotDF) = roiInds
  list(data=data,
       emmeans=plotDF,
       models=mods, method=method)
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
#' @param cft_p The cluster forming threshold or threshold for visualizing results on the p-value scale.
#' @param cft_s The cluster forming threshold or threshold for visualizing results on the RESI effect size scale.
#' @param cft_chisq The cluster forming threshold or threshold for visualizing results on the chi-squared scale.
#' @param roi Which ROI index to visualize.
#' @param index Which slices to visualize. Default is all slices
#' @param clusterID logical indicating whether to include the clusterID in the figure.
#' @param title Character passed to image.niftiImage
#' @param plane Which anatomical plane to visualize.
#' @param clusterMask logical indicating whether to mask the results with significant clusters.
#' @param alpha Adjusted p-value threshold to display clusters.
#' @param crop Whether to crop white space from image. Defaults to TRUE.
#' @param lo Whether to perform layout call to arrange images. Defaults to true.
#' @param ... Arguments passed to par.
#' @importFrom utils write.csv
#' @importFrom graphics text
#' @export
image.statMap = function(x, method=c('CEI', 'maxima', 'CMI'), cft_s=NULL, cft_p=NULL, cft_chisq=NULL, roi=NULL, index=NULL, alpha=NULL, clusterMask=TRUE, clusterID=TRUE, title='', plane=c('axial', 'sagittal', 'coronal'), crop=TRUE, lo=TRUE, ... ){
  method = tolower(method[1])
  # CFT passed as p value or effect size converted to chi-squared threshold
  if(!is.null(cft_s)){
    cft = cft_s^2 * x$sqrtSigma$n + x$sqrtSigma$df
    statmethod='S'
    thresh = cft_s
  } else if(!is.null(cft_p)){
    cft = qchisq(cft_p, df=x$sqrtSigma$df, lower.tail=FALSE)
    statmethod='p'
    thresh = -log10(cft_p)
  } else if(!is.null(cft_chisq)){
    statmethod = 'chisq'
    thresh = cft = cft_chisq
  } else { # if no cft was passed default to log10(0.05)
    cft_p = 0.05
    cft = qchisq(cft_p, df=x$sqrtSigma$df, lower.tail=FALSE)
    statmethod='p'
    thresh = -log10(cft_p)
  }
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
    image.niftiImage(stat.statMap(x, method=statmethod), BGimg = x$template, limits = thresh, index = index, plane=plane, title=title, lo=lo, crop=crop, ...)
  } else {
    if(is.character(x$mask)){
      x$mask = readNifti(x$mask)
    }
    ind = inferenceIndex(pbjObj$obsStat, method=method, cft=cft)

    st = table.statMap(x, method=method, cft_p=cft_p, cft_s=cft_s)
    imgdims = dim(stat.statMap(x, method=statmethod))
    planenum = switch(plane, "axial"=3, 'sagittal'=1, 'coronal'=2)
    otherplanes = which(!1:3 %in% planenum)
    if(!is.null(alpha)){
      # Set nonsignificant ROIs to zero
      pbjObj$ROIs[[ind]][ pbjObj$ROIs[[ind]] %in% st[which(st$`FWER p-value`>alpha),1] ] = 0
      # subset table to significant ROIs
      st = st[which(st$`FWER p-value`<=alpha),]
    } else if(!is.null(roi)){
      # set unselected ROIs to zero
      pbjObj$ROIs[[ind]][ ! pbjObj$ROIs[[ind]] %in% roi ] = 0
      # subset table to selected ROIs
      st = st[which(st[,1] %in% roi),]
    }
    # zero statMap outside of ROIs
    x$stat[pbjObj$ROIs[[ind]][ x$mask>0 ]==0] = 0

    # get coordinates from table
    coords = do.call(rbind, lapply(strsplit(st[,3], split=', '), as.numeric) )
    index = unique(coords[,planenum])
    othercoords = coords[,-planenum]

    # set graphical layout
    nCol = ceiling(sqrt(length(index)))
    nrow = ceiling(length(index)/nCol)


    par(fg='white', bg='black', mar=c(0,0,0,0), oma=c(0,0,0,0))
    # layout
    if(lo){
      lo = matrix(c(1:length(index), rep(NA, nrow * nCol - length(index) )), nrow=nrow, ncol=nCol, byrow=TRUE)
      lo[is.na(lo)] = max(lo, na.rm=TRUE)+1
      layout(lo)
    }

    if(crop){
      cr = crop(x$template, stat.statMap(x, method=statmethod))
    }

    indOrd = switch(plane[1], axial = {1:3},
                    coronal = {c(1, 3, 2)},
                    sagittal = {c(2, 3, 1)},
                    stop(paste("Orthogonal plane", plane[1], "is not valid.")))

    # if roi is non-null draw single slice with centroid of each selected ROI
    for (slic in index){
      if(clusterID){
        coordInds = which(coords[,planenum]==slic)
        coordLabels = st[coordInds,1]
        othercoords = coords[coordInds, , drop=FALSE]
        image.niftiImage(stat.statMap(x, method=statmethod), BGimg=x$template, plane=plane, index=slic, limits=thresh,
                         title=title,  lo=FALSE, crop=FALSE, ...)
        # text is offset by crop amount
        text(othercoords[,indOrd[1]]-cr$offset[indOrd[1]], othercoords[,indOrd[2]]-cr$offset[indOrd[2]], labels=coordLabels, col=par()$fg, font=2)
      } else {
        image.niftiImage(stat.statMap(x, method=statmethod), BGimg=x$template, plane=plane, index=slic-cr$offset[3], limits=thresh, title=title, lo=FALSE, crop=FALSE, ...)
      }
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
  statimg = stat.statMap(statMap, method='chisq')
  # convert Chisq to Z
  statimg[statMap$mask>0] = qnorm(pchisq(statimg[statMap$mask>0], df=statMap$sqrtSigma$df, log.p = TRUE), log.p=TRUE)
  invisible(capture.output(test <- ptfce(as.nifti(array(statimg, dim = dim(statimg))), as.nifti(array(statMap$mask, dim=dim(statMap$mask))) )))
  return(test)
}


#' Create cluster/maxima summary table
#'
#' @param x statMap object
#' @param method What statistic to provide summary for? must have run that
#' analysis using the pbjInference and mmeStat functions.
#' @param cft_s cluster forming threshold on effect size scale to show cluster p-values for.
#' @param cft_p cluster forming threshold on p-value scale to show cluster p-values for.
#' @param cft_chisq cluster forming threshold on chi-square statistic scale to show cluster p-values for.
#' @details
#' If both cft_s and cft_p are NULL, then it returns the first CFT that was used.
#'
#' @importFrom RESI chisq2S
#' @return Returns table of unadjusted and FWER adjusted p-values and other summary statistics. Results depend on what statistic
#' function was used for pbjInference.
#' @seealso [mmeStat], [cluster], [maxima], [pbjInference]
#' @export
#'
table.statMap = function(x, method=c('CEI', 'maxima', 'CMI'), cft_s=NULL, cft_p=NULL, cft_chisq=NULL){
  method = tolower(method[1])
  # CFT passed as p value or effect size converted to chi-squared threshold
  cft = cft_chisq
  if(!is.null(cft_p)){
    cft = qchisq(cft_p, df=x$sqrtSigma$df, lower.tail=FALSE)
  }
  if(!is.null(cft_s)){
    cft = cft_s^2 * x$sqrtSigma$n + x$sqrtSigma$df
  }
  statImg = stat.statMap(x, method='chisq')
  mask = x$mask
  n = x$sqrtSigma$n
  df = x$sqrtSigma$df
  x = x$pbj
  if(is.null(x)){
    if(is.null(cft) & method!='maxima'){
      stop('Must specify a cluster forming threshold\n (CFT) if pbjInference has not been run.')
    }
    x = list()
    cei = list('maxima'=FALSE, 'CEI'=FALSE, 'CMI'=FALSE)
    cei[[grep(method, tolower(names(cei))) ]] = TRUE
    x$obsStat = do.call(mmeStat, c(list(stat=statImg, mask=mask, cft=cft), cei) )
    x$ROIs = do.call(mmeStat, c(list(stat=statImg, mask=mask, cft=cft, rois=TRUE), cei) )
  }
  ind = inferenceIndex(x$obsStat, method=method, cft=cft)
  # hard-coded, but not generally applicable
  x$obsStat[[ind]] = x$obsStat[[ind]][x$obsStat[[ind]]!=0]
  if(method!='maxima'){
    Table = data.frame('Cluster ID' = 1:length(x$obsStat[[ind]]),
                       'Cluster Extent' = c(x$obsStat[[ind]]),
                       'Centroid (vox)' = sapply(1:length(x$obsStat[[ind]]), function(ind2) paste(round(colMeans(which(x$ROIs[[ind]]==ind2, arr.ind = TRUE) )), collapse=', ' )), #RNifti::voxelToWorld(
                       check.names=FALSE )
    if(method=='cmi') names(Table) = c('Cluster ID', 'Cluster Mass', 'Centroid (vox)')
    if(!is.null(x$fwerCDF)){
      if(!is.null(x$unadjCDF[[ind]])){
        Table[,'Unadjusted p-value'] =  (1-x$unadjCDF[[ind]](c(x$obsStat[[ind]])))
      }
      Table[,'FWER p-value'] = (1-x$fwerCDF[[ind]](c(x$obsStat[[ind]])))
    }
  } else {
    # 'Chi-square' assumes ROIs are indexed increasing with data frame indexing.
    Table = data.frame('Cluster ID' = 1:max(x$ROIs[[ind]]),
                       'Chi-square' = statImg[x$ROIs[[ind]]>0],
                       'Coord (vox)' = sapply(1:max(x$ROIs[[ind]]), function(ind2) paste(round(colMeans(which(x$ROIs[[ind]]==ind2, arr.ind = TRUE) )), collapse=', ' )), #RNifti::voxelToWorld(
                       check.names=FALSE )
    if(!is.null(x$fwerCDF)){
      Table[,'FWER p-value'] = (1-x$fwerCDF[[ind]](c(statImg[x$ROIs[[ind]]>0])))
    }
  }
  Table[, 'Max RESI'] = sapply(1:max(x$ROIs[[ind]]), function(ind2, statImg, ROIs, ind, n, df) chisq2S(max(statImg[ROIs[[ind]]==ind2]), df=df, n=n), statImg=statImg, ROIs=x$ROIs, ind=ind, n=n, df=df)
  # sort by the statistical value
  Table = Table[order(Table[,2], decreasing = TRUE),]
  attr(Table, 'cft') = cft
  Table
}

