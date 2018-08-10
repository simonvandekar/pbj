#' Loads NIfTI Images and Adds Synthetic Signal
#'
#' This function loads images specified in filelist and adds signal to the
#'  images. The amount and structure of signal is determined by betaimg and X.
#'  The output images are equal the images specified by filelist plus betaimg
#'  times the column of X that is not in Xred. This column is first
#'  residualized to Xred.
#' @param files a vector of .nii or .nii.gz images.
#' @param betaimg a parameter image that describes the association between X
#'  and the images in filelist. The units of this are standardized so that it
#'  is like effect size.
#' @param X the design matrix for all covariates included in the simulation
#'  model fit.
#' @param Xred the design matrix for all covariates except the column that is
#'  being multiplied by betaimg. Xred must have only one less column than X.
#' @param outfiles a vector of images to save the output.
#' @return Returns a 4d array of imaging data with synthetic signal added. The first three dimensions are equal to dim(betaimg) and the 4th dimension indexes subject.
#' @keywords power simulation
#' @importFrom abind abind
#' @importFrom stats sd
#' @export
# @examples
addSignal = function(files, betaimg, X, Xred, outfiles=NULL){

  # get column of interest
  nullinds = which(!colnames(X) %in% colnames(Xred))
  # X residualized to covariates
  x = qr.resid(qr(Xred), X[,nullinds, drop=FALSE])
  sdx = sd(x)

  # load in imaging data. Get voxelwise SD
  cat('loading images.\n')
  y = do.call(abind, list(RNifti::readNifti(files), along=4))
  cat('computing standard deviation.\n')
  sdy = apply(y, 1:3, sd)

  # load signal file
  if(is.character(betaimg)) betaimg = RNifti::readNifti(betaimg)

  # make images with signal
  y = outer(betaimg * sdy/sdx, c(x)) + y

  # write out images
  cat('writing output images.\n')
  if(is.character(outfiles[1]) & length(outfiles) == length(files))
    trash=lapply(1:length(outfiles), function(ind) RNifti::writeNifti(y[,,,ind], outfiles[ind]) )
  # return if requested
  y = y
}
