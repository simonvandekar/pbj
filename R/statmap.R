
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

statInner <- function(label, obj)
{
  if(is.null(obj))    return("")
  if(all(is.na(obj))) return(paste0(label, "NA"))

  if(class(obj)[1] == "character")  return(statFile(label, obj))
  if(class(obj)[1] == "niftiImage") return(statNifti(label, obj))
  if(class(obj)[1] == "matrix")     return(statMatrix(label, obj))

  paste0(label, "Unhandled Class(",class(obj)[1],")\n")
}

#' @export
summary.statMap <- function(object, ...)
{
  cat(paste0(
    "\nFormula: ", paste0(as.character(object$formulas[[2]]), collapse=''), paste0(as.character(object$formulas[[1]]), collapse=''), "\n",
    "\nContents:\n",
    statInner("  Stat:       ", object$stat),
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


#' Create images of a statMap
#'
#' modified from oro.nifti:::image.nifti
#'
#' @export
#' @param x the statMap object to display images of
#' @param thresh A threshold to apply to the image, defaults to 2.32
#' @param index Any selected image indexes to use for the z. defaults to NULL
#' @param col a vector of colors to use for scaled intensities defaults to a grey scale.
#' @param colpos a vector of colors to use for positive values.
#' @param colneg a vector of colors to use for negative values.
#' @param plane the plane to display, can be axial, coronal or sagittal.
#' @param xlab a title for the x axis.
#' @param ylab a title for the u axis.
#' @param axes display axes, defaults to false
#' @param oma A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text. Default to 0.
#' @param mar A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. Defaults to 0.
#' @param bg background color, defaults to black.
#' @param ... additional arguments passed to par
#' @importFrom grDevices gray
#' @importFrom graphics par
# modified from oro.nifti:::image.nifti
#' @export
image.statMap = function (x, thresh=2.32, index = NULL, col = gray(0:64/64), colpos=redyellow(64), colneg=bluecyan(64),
     plane = c("axial", "coronal", "sagittal"), xlab = "", ylab = "", axes = FALSE, oma = rep(0, 4), mar = rep(0, 4), bg = "black", ...)
  {
    object <- x
    # mask can't be empty in typical statMap object unless it's manually constructed
    if(is.null(object$mask)) object$mask=object$template
    if(is.null(object$template)) object$template=object$mask
    x = if(is.character(object$template)) readNifti(object$template) else object$template
    pixdim = RNifti::pixdim(x)
    mask = if(is.character(object$mask)) readNifti(object$mask) else object$mask
    stat = if(is.character(object$stat)) readNifti(object$stat) else object$stat

    switch(plane[1], axial = {
        aspect <- pixdim[3]/pixdim[2]
    }, coronal = {
            x <- aperm(x, c(1, 3, 2))
            mask <- aperm(mask, c(1, 3, 2))
            stat <- aperm(stat, c(1, 3, 2))
        aspect <- pixdim[4]/pixdim[2]
    }, sagittal = {
            x <- aperm(x, c(2, 3, 1))
            mask <- aperm(mask, c(2, 3, 1))
            stat <- aperm(stat, c(2, 3, 1))
        aspect <- pixdim[4]/pixdim[3]
    }, stop(paste("Orthogonal plane", plane[1], "is not valid.")))
    # permuted image dimensions

    # crop image and get image dimensions
    xinds = apply(x!=0, 1, any)
    yinds = apply(x!=0, 2, any)
    zinds = apply(x!=0, 3, any)
    x = x[xinds,,]
    x = x[,yinds,]
    x = x[,,zinds]
    stat = stat[xinds,,]
    stat = stat[,yinds,]
    stat = stat[,,zinds]
    statneg = stat
    statneg[ statneg> -thresh] = 0
    statneg = abs(statneg)
    stat[ stat<thresh ] = 0
    imgdim = dim(x)
    zlim = range(x, na.rm=TRUE)
    maxstat = max(c(stat[ stat>0], thresh), na.rm=TRUE)
    maxstatneg = max(c(statneg[ statneg>0], thresh), na.rm=TRUE)
    breaks <- c(zlim[1], seq(zlim[1], zlim[2], length = length(col) - 1), zlim[2])
    breakspos <- c(thresh, seq(thresh, maxstat, length = length(colpos)-1), maxstat)
    breaksneg <- c(thresh, seq(thresh, maxstatneg, length = length(colneg)-1), maxstatneg)
    if(is.null(index)) index = 1:imgdim[3]
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = ceiling(rep(sqrt(imgdim[3]), 2)), oma = oma, mar = mar, bg = bg)
    for (z in index) {
      # background image
      graphics::image(1:imgdim[1], 1:imgdim[2], x[, , z], col = col,
        breaks = breaks, asp = aspect, axes = axes,
        xlab = xlab, ylab = ylab, ...)
      # overlay positive
      graphics::image(1:imgdim[1], 1:imgdim[2], stat[, , z], col = colpos,
        breaks = breakspos, asp = aspect, axes = axes, add=TRUE,
        xlab = xlab, ylab = ylab, ...)
      # overlay negative
      graphics::image(1:imgdim[1], 1:imgdim[2], statneg[, , z], col = colneg,
        breaks = breaksneg, asp = aspect, axes = axes, add=TRUE,
        xlab = xlab, ylab = ylab, ...)
    }
    par(oldpar)
    invisible()
}

#' Write the statMap objects out
#'
#' Given a statMap object and a directory write the objects as stat.nii.gz, sqrtSigma.nii.gz and summary.txt
#' @param x the statMap object to write out
#' @param outdir the directory to write into
#' @param ... additional arugment, unused
#' @return a list of what was written
#' @export
write.statMap <- function(x,outdir, ...)
{
  statimg  = file.path(outdir, 'stat.nii.gz')
  resimg   = file.path(outdir, 'sqrtSigma.nii.gz')
  summaryf = file.path(outdir, 'summary.txt')
  if(is.character(x$stat)){
    file.copy(x$stat, statimg)
    file.copy(x$sqrtSigma, resimg)
  } else {
    cat('Writing output images.\n')
    dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
    writeNifti(x$stat, statimg)

    cat('Writing sqrtSigma 4d image.\n')
    # reshape res matrix into 4d image
    # memory intensive
    # overwriting res
    x$sqrtSigma = lapply(1:ncol(x$sqrtSigma), function(ind){ x$mask[ x$mask==1] = x$sqrtSigma[,ind]; x$mask} )
    # combine into 4d array
    x$sqrtSigma = do.call(abind::abind, c(x$sqrtSigma, along=4))
    writeNifti(updateNifti(x$sqrtSigma, x$mask), resimg)
    res = resimg
  }
  return(list(stat=statimg, sqrtSigma=resimg, summary=summaryf))
}

