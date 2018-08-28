
sizing <- function(fn) paste0(round(file.info(fn)$size/1024^2, 2), "M")
statFile  <- function(label, fn) paste0(label, fn, "' ", sizing(fn), "\n")
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
  
  paste0(label, "Unhandled Class(",class(obj)[1],")")
}

#' @export
summary.statMap <- function(object, ...)
{
  cat(paste0(
    "\nFormula: ", paste0(as.character(object$formulas[[2]]), collapse=''), paste0(as.character(object$formulas[[1]]), collapse=''), "\n",
    "\nContents:\n",
    statInner("  Stat:       '", object$stat),
    statInner("  Sqrt Sigma: '", object$sqrtSigma),
    statInner("  Mask:       '", object$mask),
    statInner("  Template:   '", object$template),
    "  Robust:     ", object$robust, "\n"
  ))
}

#' @export
print.statMap <- function(x, ...)
{
  cat(summary(x, ...))
}

plot.statMap <- function(x, slice=1, ...)
{
  stop("FILL IN PLOTTING FUNCTION HERE. SLICE IS an example of an additional parameter")
}

redyellow = colorRampPalette(c('red', 'yellow'))
bluecyan = colorRampPalette(c('blue', 'cyan'))

# modified from oro.nifti:::image.nifti
image.statMap = function (statmap, thresh=2.32, index = NULL, col = gray(0:64/64), colpos=redyellow(64), colneg=bluecyan(64),
     plane = c("axial", "coronal", "sagittal"), xlab = "", ylab = "", axes = FALSE, oma = rep(0, 4), mar = rep(0, 4), bg = "black", ...) 
{
    x = if(is.character(statmap$template)) readNifti(statmap$template) else statmap$template
    pixdim = RNifti::pixdim(x)
    mask = if(is.character(statmap$mask)) readNifti(statmap$mask) else statmap$mask
    stat = if(is.character(statmap$stat)) readNifti(statmap$stat) else statmap$stat

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
    xinds = apply(mask==1, 1, any)
    yinds = apply(mask==1, 2, any)
    zinds = apply(mask==1, 3, any)
    x = x[xinds,,]
    x = x[,yinds,]
    x = x[,,zinds]
    stat = stat[xinds,,]
    stat = stat[,yinds,]
    stat = stat[,,zinds]
    statneg = stat
    statneg[ statneg>= -thresh] = 0
    statneg = abs(statneg)
    stat[ stat<=thresh ] = 0 
    imgdim = dim(x)
    zlim = range(x, na.rm=TRUE)
    maxstat = max(c(stat, thresh), na.rm=TRUE)
    maxstatneg = max(c(statneg, thresh), na.rm=TRUE)
    breaks <- c(zlim[1], seq(zlim[1], zlim[2], length = length(col) - 1), zlim[2])
    breakspos <- c(thresh, seq(thresh, maxstat, length = length(colpos) - 1), maxstat)
    breaksneg <- c(thresh, seq(thresh, maxstatneg, length = length(colneg) - 1), maxstatneg)
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
