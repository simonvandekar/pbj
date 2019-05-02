#' Create images of a NiFTI object
#'
#' modified from oro.nifti:::image.nifti
#'
#' @export
#' @param x the Nifti object to display images of
#' @param bgimg background image to use.
#' @param thresh A threshold to apply to the image, defaults to 0
#' @param index Any selected image planes. defaults to NULL
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
image.niftiImage = function (x, bgimg=NULL, thresh=0, index = NULL, col = gray(0:64/64), colpos=redyellow(64), colneg=bluecyan(64),
                          plane = c("axial", "coronal", "sagittal"), xlab = "", ylab = "", axes = FALSE, oma = rep(0, 4), mar = rep(0, 4), bg = "black", ...)
{
  eps = 10^-6
  thresh = thresh + eps
  object <- x
  stat = if(is.character(object)) readNifti(object)  else object
  x = if(is.character(bgimg)){
      readNifti(bgimg)
    } else if(is.null(bgimg)){
      bgimg = stat
      } else bgimg
  pixdim = RNifti::pixdim(x)

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
  par(mfrow = ceiling(rep(sqrt(length(index)), 2)), oma = oma, mar = mar, bg = bg)
  for (z in index) {
    # background image
    graphics::image(1:imgdim[1], 1:imgdim[2], x[, , z], col = col,
                    breaks = breaks, asp = aspect, axes = axes,
                    xlab = xlab, ylab = ylab, ...)
    # overlay positive
    if(thresh!=maxstat)
    graphics::image(1:imgdim[1], 1:imgdim[2], stat[, , z], col = colpos,
                    breaks = breakspos, asp = aspect, axes = axes, add=TRUE,
                    xlab = xlab, ylab = ylab, ...)
    # overlay negative
    if(thresh != maxstatneg)
    graphics::image(1:imgdim[1], 1:imgdim[2], statneg[, , z], col = colneg,
                    breaks = breaksneg, asp = aspect, axes = axes, add=TRUE,
                    xlab = xlab, ylab = ylab, ...)
  }
  par(oldpar)
  invisible()
}
