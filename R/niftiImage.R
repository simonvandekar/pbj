#' Create images of a NiFTI object
#'
#' modified from oro.nifti:::image.nifti
#'
#' @export
#' @param x the Nifti object to display images of
#' @param bgimg background image to use.
#' @param thresh A lower (an optionally upper length 2 vector) threshold to apply to the image, defaults to 0
#' @param index Any selected image planes. defaults to NULL
#' @param col a vector of colors to use for scaled intensities defaults to a grey scale.
#' @param colpos a vector of colors to use for positive values.
#' @param colneg a vector of colors to use for negative values.
#' @param plane the plane to display, can be axial, coronal or sagittal.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param axes display axes, defaults to false
#' @param oma A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text. Default to 0.
#' @param mar A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. Defaults to 0.
#' @param bg background color, defaults to black.
#' @param other A function used to add features to the plot.
#' @param ... additional arguments passed to par
#' @importFrom grDevices gray
#' @importFrom graphics par
#' @importFrom RNifti readNifti
#' @export
#' @examples
#' library(RNifti)
#' library(pain21)
#' pain = pain21()
#' templ = readNifti(pain$template)
#' # image(templ)
#' # image(templ, plane='coronal')
#'
image.niftiImage = function (x, bgimg = NULL, thresh = 0, index = NULL, col = gray(0:64/64),
                             colpos = pbj:::redyellow(64), colneg = pbj:::bluecyan(64), plane = c("axial",
                                                                                                  "coronal", "sagittal"), xlab = "", ylab = "", axes = FALSE,
                             oma = rep(0, 4), mar = rep(0, 4), bg = "black", other=function(){}, ...)
{
  eps = 10^-6
  thresh = thresh + eps
  object <- x
  stat = if (is.character(object))
    readNifti(object)
  else object
  if (is.character(bgimg)) {
    x=readNifti(bgimg)
  }
  else if (is.null(bgimg)) {
    x = stat
    # if only one argument was passed, display as the background img (foreground is blank)
    if(thresh==eps){
     thresh=max(stat)+eps
    }
  }
  else x=bgimg
  pixdim = RNifti::pixdim(x)
  switch(plane[1], axial = {
    aspect <- pixdim[3]/pixdim[2]
  }, coronal = {
    x <- aperm(x, c(1, 3, 2))
    #mask <- aperm(mask, c(1, 3, 2))
    stat <- aperm(stat, c(1, 3, 2))
    aspect <- pixdim[4]/pixdim[2]
  }, sagittal = {
    x <- aperm(x, c(2, 3, 1))
    #mask <- aperm(mask, c(2, 3, 1))
    stat <- aperm(stat, c(2, 3, 1))
    aspect <- pixdim[4]/pixdim[3]
  }, stop(paste("Orthogonal plane", plane[1], "is not valid.")))
  xinds = apply(x != 0, 1, any)
  yinds = apply(x != 0, 2, any)
  zinds = apply(x != 0, 3, any)
  x = x[xinds, , ]
  x = x[, yinds, ]
  x = x[, , zinds]
  stat = stat[xinds, , ]
  stat = stat[, yinds, ]
  stat = stat[, , zinds]
  statneg = stat
  statneg[statneg > -thresh[1]] = 0
  statneg = abs(statneg)
  stat[stat < thresh[1]] = 0
  imgdim = dim(x)
  zlim = range(x, na.rm = TRUE)

  maxstat = max(replace(c(stat[stat > 0], thresh), !is.finite(c(stat[stat >
                                                                       0], thresh)), NA), na.rm = TRUE)
  maxstatneg = max(replace(c(statneg[statneg > 0], thresh),
                           !is.finite(c(statneg[statneg > 0], thresh)), NA), na.rm = TRUE)

  breaks <- c(zlim[1], seq(zlim[1], zlim[2], length = length(col) -
                             1), zlim[2])
  if(length(thresh)<2){
    breakspos <- c(thresh[1], seq(thresh[1], maxstat, length = length(colpos) -
                                    1), maxstat)
    breaksneg <- c(thresh[1], seq(thresh[1], maxstatneg, length = length(colneg) -
                                    1), maxstatneg)
  } else {
    breakspos <- c(thresh[1], seq(thresh[1], thresh[2], length = length(colpos) -
                                    1), maxstat)
    breaksneg <- c(thresh[1], seq(thresh[1], thresh[2], length = length(colneg) -
                                    1), maxstatneg)
  }
  if (is.null(index)) {
    index = 1:imgdim[3]
    par(mfrow = ceiling(rep(sqrt(length(index)), 2)), oma = oma,
        mar = mar, bg = bg)
    oldpar <- par(no.readonly = TRUE)
  }
  for (z in index) {
    graphics::image(1:imgdim[1], 1:imgdim[2], x[, , z], col = col,
                    breaks = breaks, asp = aspect, axes = axes, xlab = xlab,
                    ylab = ylab, ...)
    if (thresh[1] != maxstat)
      graphics::image(1:imgdim[1], 1:imgdim[2], stat[,
                                                     , z], col = colpos, breaks = breakspos, asp = aspect,
                      axes = axes, add = TRUE, xlab = xlab, ylab = ylab,
                      ...)
    if (thresh[1] != maxstatneg)
      graphics::image(1:imgdim[1], 1:imgdim[2], statneg[,
                                                        , z], col = colneg, breaks = breaksneg, asp = aspect,
                      axes = axes, add = TRUE, xlab = xlab, ylab = ylab,
                      ...)
  }
  other()
  # if (length(index) != 1)
  #   par(oldpar)
  invisible()
}
