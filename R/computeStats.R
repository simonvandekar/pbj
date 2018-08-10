#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
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
#' @param Winv Inverse weights for regression model. Inverse of W. Required when
#' passing variance images as the inverse weights.
#' @param robust Compute robust standard error estimates? Defaults to TRUE.
#'   Uses HC3 SE estimates from REF.
#' @param statfile nii or nii.gz file to save out 3d statistical image.
#' @param resfile nii or nii.gz file to save out 4d covariance image.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @importFrom abind abind
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical nifti object. If ncol(X) = ncol(Xred)+1, then this is a Z-statistic map, otherwise it is a chi^2-statistic map.}
#'   \item{res}{The 4d covariance object. This is a V by n matrix R, such that R \%*\% t(R) = Sigma.}
#' }
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @export
computeStats = function(files, X, Xred, Xfiles=NULL, mask, W=rep(1, nrow(X)), Winv=NULL, robust=TRUE, statfile=NULL, resfile=NULL, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001

  # check if inverse weights are given
  if(!is.null(Winv)){
    W = Winv
    Winv = TRUE
  } else {
    Winv = FALSE
  }

  if(is.character(files)){
    n=length(files)
    if(nrow(X)!=n)
      stop('length(files) and nrow(X) must be the same.')
    res = do.call(abind::abind, list(RNifti::readNifti(files), along=4))
    } else {
      n = nrow(X)
      res = files
      rm(files)
    }

  # load mask
  if(is.character(mask))
    mask = RNifti::readNifti(mask)

  # load images
  res = t(apply(res, 4, function(x) x[mask==1]))

  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)
  if(df>1 & robust)
    stop('Robust covariance only available for testing a single parameter.')

  if(is.character(W)){
    cat('Weights are voxel-wise.\n')
    voxwts = TRUE
    W = do.call(abind::abind, list(RNifti::readNifti(W), along=4))
    W = t(apply(W, 4, function(x) x[mask==1]))
    W = sqrt(W)
  } else {
    voxwts=FALSE
    # compute W half
    W = c(sqrt(W))
  }
  # inverse weights were passed
  if(Winv) W[W!=0] = 1/W[W!=0]
  X1 = X[,peind]
  # this is a pointwise matrix multiplication if W was passed as images
  res = res * W

  # fit model to all image data
  # if weights are voxel specific then design must also be treated separately
  if(voxwts){
    cat('Running voxel-wise weighted linear models.\n')

    if(!robust){
      if(df==1){
        # cannot be parallelized due to memory use
        qrs = apply(W, 2, function(Wcol) qr(X * Wcol))
        # This should work in parallel
        seX1 = sqrt(do.call(c, parallel::mclapply(qrs, function(qrval) chol2inv(qr.R(qrval))[peind,peind], mc.cores=mc.cores)))
        # cannot be parallelized due to memory use
        num = do.call(c, lapply(1:ncol(res), function(ind) qr.coef(qrs[[ind]], res[,ind])[peind]) )
        # cannot be parallelized due to memory use
        res = do.call(rbind, lapply(1:ncol(res), function(ind) qr.resid(qrs[[ind]], res[,ind])) )
        rm(qrs, W)
      } else {
        num = do.call(c, lapply(1:ncol(res), function(ind) sum(qr.resid(qr(Xred * W[,ind]), res[,ind])^2)) )
        res = do.call(rbind, lapply(1:ncol(res), function(ind) qr.resid(qr(X * W[,ind]), res[,ind])) )
      }
    }

    if(robust){
      res = lapply(1:ncol(res), function(ind) lm(res[,ind] ~ -1 + I(X * W[,ind]), model=FALSE) )

      # get parameter estimates
      stat = do.call(rbind, parallel::mclapply(res, coefficients, mc.cores=mc.cores ))[,peind]

      cat('Getting voxel-wise hat values.\n')
      h = do.call(rbind, parallel::mclapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}, mc.cores=mc.cores ))

      cat('Getting voxel-wise residuals for covariate and outcome vectors.\n')

      res = do.call(rbind, parallel::mclapply(res, residuals, mc.cores=mc.cores))
      X1res = do.call(rbind, lapply(1:nrow(res), function(ind) qr.resid(qr(Xred * W[,ind]), X1 * W[,ind])) )
      res = res * X1res /(1-h)
      A = rowSums(X1res^2)
      rm(h, X1res)
    }

  # else weights are the same for all voxels
  } else {
    if(!robust){
      if(df==1){
	QR = qr(X * W)
        num = qr.coef(QR, res)[peind,]
        seX1 = sqrt(chol2inv(qr.R(QR))[peind,peind])
        res = t(qr.resid(QR, res))
	rm(QR)
      } else {
        num = colSums(qr.resid(qr(Xred * W), res)^2)
        res = t(qr.resid(qr(X * W), res))
      }
    }

    if(robust){
      res = lm(res ~ -1 + I(X * W), model=FALSE)

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
    # assume t-statistics if df==1
    if(df==1){
      stattemp = sqrt(stattemp/rdf)
      stattemp = num/stattemp / seX1
      # convert to z-statistics
      stattemp = qnorm(pt(stattemp, df=rdf))
    } else {
      num = num - stattemp
      stattemp = num/stattemp * rdf/df
      # convert to chisquared
      stattemp = qchisq(pf(stattemp, df1=df, df2=rdf), df=df)
    }
    stat = mask
    stat[ mask==1] = stattemp
  }

  if(robust){
    stattemp = stat*A/sqrt(rowSums(res^2))
    stat = mask
    stat[ stat==1] = stattemp
  }

  if(!is.null(resfile)){
    # need to re-form residuals into an array first
    writeNifti(res, resfile)
  }

  if(!is.null(statfile))
    writeNifti(stat, statfile)

  # returns if requested
  out = list(stat=stat, res=res)
}
