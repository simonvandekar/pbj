#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and covariance matrix which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param files Character vector of subject images to be modeled as an outcome
#'  variable OR 4d array of imaging data OR 4d nifti object.
#' @param X Design matrix for full model.
#' @param Xred Design matrix for reduced model. If robust=TRUE then this must
#'  have one less column than X.
#' @param Xfiles n X p matrix of images where n is the number of subjects and
#'  each column corresponds to an imaging covariate. Currently, not available.
#' @param mask File name for a binary mask file.
#' @param W Weights for regression model. Can be used to deweight noisy
#'  observations. Same as what should be passed to lm.
#' @param robust Compute robust standard error estimates? Defaults to TRUE.
#'   Uses HC3 SE estimates from REF.
#' @param statfile nii or nii.gz file to save out 3d statistical image.
#' @param resfile nii or nii.gz file to save residuals.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @export
#' @examples
computeStats = function(files=NULL, X=NULL, Xred=NULL, Xfiles=NULL, mask=NULL, W=rep(1, nrow(X)), robust=TRUE, statfile=NULL, resfile=NULL, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001
  if(any( is.null(list(files, X, Xred, mask ) )))
    stop('One or more required arguments unspecified.')

  if(is.character(files)){
    n=length(files)
    if(nrow(X)!=n)
      stop('length(files) and nrow(X) must be the same.')
    } else {
      res = files
      rm(files)
    }

  # load mask
  if(is.character(mask))
    mask = RNifti::readNifti(mask)

  # load images
  res = do.call(abind, list(RNifti::readNifti(files), along=4))
  res = t(apply(res, 4, function(x) x[mask==1]))

  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)
  if(df>1 & robust)
    stop('Robust covariance only available for testing a single parameter.')

  if(is.character(W)){
    cat('Weights are voxel-wise.\n')
    voxwts = TRUE
    W = do.call(abind, list(RNifti::readNifti(W), along=4))
    W = apply(W, 4, function(x) x[mask==1])
    W = sqrt(W)
  } else {
    voxwts=FALSE
    # compute W half
    W = c(sqrt(W))
    X = X
  }
  X1 = X[,peind]
  # this is a pointwise matrix multiplication if W was passed as images
  res = res * W

  # fit model to all image data
  # if weights are voxel specific then design must also be treated separately
  if(voxwts){
    cat('Running voxel-wise weighted linear models.\n')

    if(!robust){
      num = do.call(c, mclapply(1:ncol(res), function(ind) sum(qr.resid(qr(Xred * W[,ind]), res[,ind])^2), mc.cores=mc.cores) )
      res = do.call(rbind, mclapply(1:ncol(res), function(ind) qr.resid(qr(X * W[,ind]), res[,ind]), mc.cores=mc.cores))
    }

    if(robust){
      res = mclapply(1:ncol(res), function(ind) lm(res[,ind] ~ -1 + I(X * W[,ind]), model=FALSE), mc.cores=mc.cores)

      # get parameter estimates
      stat = do.call(rbind, mclapply(res, coefficients, mc.cores=mc.cores ))[,peind]

      cat('Getting voxel-wise hat values.\n')
      h = do.call(rbind, mclapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}, mc.cores=mc.cores ))

      cat('Getting voxel-wise residuals for covariate and outcome vectors.\n')
      res = do.call(rbind, mclapply(res, residuals, mc.cores=mc.cores))
      X1res = do.call(rbind, mclapply(1:ncol(res), function(ind) qr.resid(qr(Xred * W[,ind]), X1 * W[,ind]), mc.cores=mc.cores ))
      res = t(res * X1res /(1-h))
      A = rowSums(X1res^2)
      rm(h, X1res)
    }

  # else weights are the same for all voxels
  } else {
    if(!robust){
      num = colSums(qr.resid(qr(Xred * W), res)^2)
      res = t(qr.resid(qr(X * W), res))
    }

    if(robust){
      res = lm(res ~ -1 + (X * W), model=FALSE)

      # get parameter estimates
      stat=coefficients(res)[peind,]

      # compute hat values
      h = rowSums(qr.Q(res$qr)^2)
      h = ifelse(h>=1, 1-eps, h)

      # get residuals
      res = residuals(res)

      # residualize variable of interest to covariates
      X1res = qr.resid(qr(Xred * W), X1 * W) # Formula XX in the paper
      A = sum(X1res^2)
      # compute half of covariance of parameter of interest
      # divides by 1-h to use the HC3 version discussed by Long and Ervin
      # https://pdfs.semanticscholar.org/1526/72b624b44b12250363eee602554fe49ca782.pdf
      res = t(res *  (X1res/(1-h)))
    }
  } # end if-else voxwts

  # compute statistical image
  if(!robust){
    stattemp = rowSums(res^2)
    num = num - stattemp
    stattemp = num/stattemp * rdf/df
    stat = mask
    stat[ mask==1] = stattemp
  }

  if(robust){
    stattemp = stat*A/sqrt(rowSums(res^2))
    stat = mask
    stat[ stat==1] = stattemp
  } else {
    # compute non-robust statistical image
  }

  if(!is.null(resfile)){
    # need to reform residuals into an array first
    writeNifti(res, resfile)
  }

  if(!is.null(statfile))
    writeNifti(stat, statfile)

  # returns if requested
  out = list(stat=stat, res=res)
}
