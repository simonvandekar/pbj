#' Prepares Simulation Data for Bootstrapping
#'
#' Residualizes the images in files to the model form and writes the output to
#' outfiles optionally smoothes with sm (in mm FWHM) using susan prior to
#' residualizing if smoutfiles is specified smoothed images are saved into that
#' directory. Used to run simulations to assess power and type 1 error for
#' papers. This creates images that are residualized to the covariates which
#' can then be bootstrapped to generate a sample where there is the potential
#' for heteroskedasticity/nonexchangeability, but where the covariates are
#' unassociated with the mean of the outcome.
#' @param files Character vector of subject images to be modeled as an outcome
#'  variable.
#' @param form mgcv or lm style formula.
#' @param dat Data frame containing covariates used by form.
#' @param mask Character giving location of mask image.
#' @param outfiles Character vector of images to save the residual output images.
#' @param sm Numeric giving the smoothing amount in mm FWHM to perform before creating residuals. Smoothing is performed using fsl's susan.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @return No returned value. This functions saves out nifti images files after residualizing to the model specified by form and dat. The residuals of files are saved as the corresponding element in outfiles.
#' @keywords power simulation, parametric bootstrap, type 1 error simulations, null simulations
#' @importFrom abind abind
#' @export
# @examples
simulationSetup = function(files=NULL, form=NULL, dat=NULL, mask=NULL, outfiles=NULL, smoutfiles=NULL, sm=0){
  if(any( is.null(list(files, form, dat, mask, outfiles ) )))
    stop('One or more required arguments unspecified.')

  # smooth using susan
  if(sm>0){
    smsigma =  sm/sqrt(8 * log(2))
    cat('smoothing files with susan\n')
    # writes output if outfiles given
    if(is.character(smoutfiles[1])){
      dir.create(dirname(smoutfiles[1]), showWarnings=FALSE)
      y = mclapply(1:length(files), function(ind) fslr::susan(files[ind], smoutfiles[ind], retimg=TRUE, sigma=smsigma) )
      # else returns 4d array
    } else {
      y = mclapply(1:length(files), function(ind) fslr::susan(files[ind], retimg=TRUE, sigma=smsigma) )
      y = do.call(abind::abind, list(y, along=4))
    }
  } else {
    cat('loading images.\n')
    y = do.call(abind::abind, list(readNifti(files), along=4))
  }

    # run linear model to get residuals
    if(is.character(mask)){
      mask = RNifti::readNifti(mask)
    }
    y = t(apply(y, 4, function(x) x[mask==1]))
    X = getDesign(form, data=dat)
    cat('regressing out covariates.\n')
    y = qr.resid(qr(X), y)

  # assumes your putting all files in the same directory or that all directories exist
  dir.create(dirname(outfiles[1]), showWarnings=FALSE)

  # save out images
  temp = mask
  trash = lapply(1:nrow(y), function(ind){ temp[ temp==1] = y[ind,]
                        RNifti::writeNifti(temp, outfiles[ind])
                        })
  return(NULL)
}
