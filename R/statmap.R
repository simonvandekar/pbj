
sizing <- function(fn) paste0(round(file.info(fn)$size/1024^2, 2), "M")
statFile  <- function(label, fn) paste0(label, "'", fn, "' ", sizing(fn), "\n")
statNifti <- function(label, im)
{
  paste0(
    label, "nifti[",
    paste0(dim(im), collapse=" x "),
    "] Pixel dimensions: ",
    paste0(attr(im, "pixdim"), c(attr(im, "pixunits")[1], attr(im, "pixunits")), collapse=" x "),
    '\n'
  )
}
statMatrix <- function(label, mm)
{
  paste0(
    label, "matrix[",
    paste0(dim(mm), collapse=" x "),
    "]\n"
  )
}

statVector <- function(label, mm)
{
  paste0(
    label, "vector[",
    length(mm),
    "]\n"
  )
}

statInner <- function(label, obj)
{
  if(is.null(obj))    return("")
  if(all(is.na(obj))) return(paste0(label, "NA"))

  if(class(obj)[1] == "character")  return(statFile(label, obj))
  if(class(obj)[1] == "niftiImage") return(statNifti(label, obj))
  if(class(obj)[1] == "matrix")     return(statMatrix(label, obj))
  if(class(obj)[1] == "matrix")     return(statVector(label, obj))

  paste0(label, "Unhandled Class(",class(obj)[1],")\n")
}

#' @export
summary.statMap <- function(object, ...)
{
  cat(paste0(
    "\nFormula: ", paste0(as.character(object$formulas[[2]]), collapse=''), paste0(as.character(object$formulas[[1]]), collapse=''), "\n",
    "\nContents:\n",
    statInner("  Stat:       ", object$stat),
    statInner("  Coef:       ", object$coef),
    statInner("  Sqrt Sigma: ", object$sqrtSigma),
    statInner("  Mask:       ", object$mask),
    statInner("  Template:   ", object$template),
    "  Robust:     ", object$robust, "\n"
  ))
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
#' @return a list of what was written
#' @export
write.statMap <- function(x, outdir){
  statimg  = file.path(outdir, 'stat.nii.gz')
  coefimg   = file.path(outdir, 'coef.nii.gz')
  res   = file.path(outdir, 'sqrtSigma.rds')
  if(is.character(x$stat)){
    file.copy(x$stat, statimg)
    file.copy(x$sqrtSigma, res)
    file.copy(x$coef, coefimg)
  } else {
    message('Writing stat and coef images.')
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    writeNifti(stat.statMap(x), statimg)
    writeNifti(coef.statMap(x), coefimg)

    message('Writing sqrtSigma object.')
    saveRDS(x$sqrtSigma, file = res)
  }


  pbj = x$pbj
  if(!is.null(pbj)){
    #if(is.null(x$template)) x$template = x$mask
    #sform = do.call(rbind, RNifti::niftiHeader(x$template)[c('srow_x', 'srow_y', 'srow_z')])
    #voxvol = prod(RNifti::pixdim(x$template))
    inference=list()
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
  coef = simplify2array(lapply(1:nrow(object$coef), function(coefv){ object$mask[object$mask!=0] = coefv; object$mask}))
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
