#' Create images of a NiFTI object
#'
#' modified from oro.nifti:::image.nifti
#'
#' @export
#' @param x the Nifti object to display images of
#' @param bgimg background image to use.
#' @param limits A lower (and optionally upper length 2 vector) limitsold to apply to the image, defaults to 0
#' @param nrow number of rows to display the slices as. NULL is default and it's a square shape determined by length of `index` argument.
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
#' @param fg foreground color, defaults to white.
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
image.niftiImage = function (x, bgimg = NULL, limits = 0, nrow=NULL, index = NULL, col = gray(0:64/64),
                             colpos = pbj:::redyellow(64), colneg = pbj:::bluecyan(64), plane = c("axial",
                                                                                                  "coronal", "sagittal"),
                             xlab = "", ylab = "", axes = FALSE,
                             oma = rep(0, 4), mar = rep(0, 4), bg = "black",  fg='white', other=function(){}, ...)
{
  eps = 10^-6
  limits = limits + eps
  object <- x
  stat = if (is.character(object))
    readNifti(object)
  else object
  if (is.character(bgimg)) {
    x=readNifti(bgimg)
  }
  else if (is.null(bgimg)) {
    # this sets foreground image to be the background too
    x = stat
    # if only one argument was passed, display as the background img (foreground is blank)
    # this sets the limitsold high enough so that no overlay is shown (if it wasn't set already)
    if(limits[1]==eps){
     limits[1]=max(abs(stat))+eps
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
  statneg[statneg > -limits[1]] = 0
  statneg = abs(statneg)
  stat[stat < limits[1]] = 0
  imgdim = dim(x)
  zlim = range(x, na.rm = TRUE)

  maxstat = max(replace(c(stat[stat > 0], limits), !is.finite(c(stat[stat >
                                                                       0], limits)), NA), na.rm = TRUE)
  maxstatneg = max(replace(c(statneg[statneg > 0], limits),
                           !is.finite(c(statneg[statneg > 0], limits)), NA), na.rm = TRUE)

  breaks <- c(zlim[1], seq(zlim[1], zlim[2], length = length(col) -
                             1), zlim[2])
  if(length(limits)<2){
    breakspos <- c(limits[1], seq(limits[1], maxstat, length = length(colpos) -
                                    1), maxstat)
    breaksneg <- c(limits[1], seq(limits[1], maxstatneg, length = length(colneg) -
                                    1), maxstatneg)
  } else {
    breakspos <- c(limits[1], seq(limits[1], limits[2], length = length(colpos) -
                                    1), maxstat)
    breaksneg <- c(limits[1], seq(limits[1], limits[2], length = length(colneg) -
                                    1), maxstatneg)
  }
  if (is.null(index)) {
    index = 1:imgdim[3]
  }
  if(is.null(nrow)){
    nCol = ceiling(sqrt(length(index)))
    nrow = ceiling(length(index)/nCol)
  } else {
    nCol = ceiling(length(index)/nrow)
  }
  # layout
  lo = matrix(c(1:length(index), rep(NA, nrow * nCol - length(index) )), nrow=nrow, ncol=nCol, byrow=TRUE)
  # if(colorbar){
  #   lo = cbind(lo, max(lo, na.rm=TRUE)+1)
  #   lo[is.na(lo)] = max(lo, na.rm=TRUE)+1
  #   layout(lo, widths=c(rep(1, nrow), 0.125))
  # } else {
    lo[is.na(lo)] = max(lo, na.rm=TRUE)+1
    layout(lo)
  # }

  par(oma = oma,
      mar = mar, bg = bg)
  oldpar <- par(no.readonly = TRUE)
 for(z in index){
      graphics::image(1:imgdim[1], 1:imgdim[2], x[, , z], col = col,
                      breaks = breaks, asp = aspect, axes = axes, xlab = xlab,
                      ylab = ylab, ...)
      if (limits[1] != maxstat)
        graphics::image(1:imgdim[1], 1:imgdim[2], stat[,
                                                       , z], col = colpos, breaks = breakspos, asp = aspect,
                        axes = axes, add = TRUE, xlab = xlab, ylab = ylab,
                        ...)
      if (limits[1] != maxstatneg)
        graphics::image(1:imgdim[1], 1:imgdim[2], statneg[,
                                                          , z], col = colneg, breaks = breaksneg, asp = aspect,
                        axes = axes, add = TRUE, xlab = xlab, ylab = ylab,
                        ...)
  }
  # test = plot_grid(plotlist=res, nrow=nrow, ncol=nCol)
  # cb = function(){
  #   par(oma=c(0,0,0,0), mar=c(8, 4, 8, 0.5), mgp=c(3,0.6,0), fg=fg, col.axis=fg, col.lab=fg, col.main = fg, col.sub=fg, bg=bg)
  #   colorBar(pbj:::redyellow(64), min=limits, max=maxstat, nticks=4, ylab = barLab)
  # }
  # plot_grid(test, cb, ncol = 2, rel_widths = c(6,1))
  # browser()
  # # height of plot device
  # wh = dev.size('in')
  # # height of one brain image
  # barHeight = wh[2]/nrow
  # marginHeight = (wh[2]-barHeight)/2
  # # relative width of column for colorBar
  # cbWidth = wh[1] * 0.125/(nrow+0.125)
  # #par(oma=c(0,0,0,0), mai=c(marginHeight, cbWidth*0.1, marginHeight, cbWidth*0.02), mgp=c(3,0.6,0), fg=fg, col.axis=fg, col.lab=fg, col.main = fg, col.sub=fg)
  # colorBar(pbj:::redyellow(64), min=limits, max=maxstat, nticks=4, ylab = barLab)
  # other()
  # invisible()
}
