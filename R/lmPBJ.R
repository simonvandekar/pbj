#' Computes Statistical Map for Neuroimaging Data
#'
#' This function computes a statistical map and residuals which are the
#'  objects necessary to perform the parametric bootstrap joint (PBJ) inference
#'  procedure.
#' @param images Character vector of subject images to be modeled as an outcome
#'  variable OR 4d array of imaging data OR 4d nifti object.
#' @param form formula or character that can be coerced to a formula or a design
#' matrix for the full model.
#' @param formred formula or character that can be coerced to a formula or a
#' design matrix for reduced model. If robust=TRUE then this must have one less
#' column than X.
#' @param mask File name for a binary mask file or niftiImage object.
#' @param data R data frame containing variables in form. If form and formred
#' are matrices then this can be NULL.
#' @param W Weights for regression model. Can be used to deweight noisy
#'  observations. Same as what should be passed to lm.
#' @param Winv Inverse weights for regression model. Inverse of W. Required when
#' passing variance images as the inverse weights.
#' @param template Template image used for visualization.
#' @param formImages n X p matrix of images where n is the number of subjects and
#'  each column corresponds to an imaging covariate. Currently, not supported.
#' @param robust Logical, compute robust standard error estimates? Defaults to TRUE.
#'   Uses HC3 SE estimates from Long and Ervin 2000.
#' @param sqrtSigma Logical: return V X n matrix sqrtSigma? Defaults to TRUE (described below).
#' Required to use pbj function.
#' @param transform Logical: use transformation from equation (5) of Vandekar et al. 2019 (Biometrics)?
#' Defaults to TRUE.
#' (instead of niftiImage objects; defaults to FALSE).
#' @param outdir If specified, output is saved as NIfTIs and statMap object is
#' saved as strings. This approach conserves memory, but has longer IO time.
#' Currently, not supported.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @importFrom abind abind
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical values where mask!=0. If ncol(X) = ncol(Xred)+1, then this is a Z-statistic map, otherwise it is a chi^2-statistic map.}
#'   \item{coef}{A 4d niftiImage giving the parameter estimates for covariates only in the full model.}
#'   \item{sqrtSigma}{The covariance object. This is a V by n matrix R, such that R \%*\% t(R) = hatSigma.}
#'   \item{mask}{The input mask.}
#'   \item{template}{The background template used for visualization.}
#'   \item{formulas}{A list containing the full and reduced models.}
#'   \item{robust}{A logical indicating the input setting.}
#'   \item{df}{The numerator degrees of freedom of the test.}
#'   \item{rdf}{The numerator degrees of freedom of the test.}
#' }
#' stat=stat, sqrtSigma=res, mask=mask, template=template, formulas=list(form, formred), robust=robust, df=df, rdf=rdf
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @export
lmPBJ = function(images, form, formred, mask, data=NULL, W=NULL, Winv=NULL, template=NULL, formImages=NULL, robust=TRUE, sqrtSigma=TRUE, transform=TRUE, outdir=NULL, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001

  if(!is.matrix(form) & !is.matrix(formred)){
    X = getDesign(form, data)
    Xred = if(!is.null(formred)) getDesign(formred, data) else NULL
  } else {
    X = form
    Xred = formred
    form <- formred <- NULL

  }

  if(class(images)[1] != 'niftiImage'){
    n=length(images)
    images = as.character(images)
    images = gsub(" +$", "", images)
    if(nrow(X)!=n)
      stop('length(images) and nrow(X) must be the same.')
    res = do.call(abind::abind, list(RNifti::readNifti(images), along=4))
    } else {
      n = nrow(X)
      res = images
      rm(images)
    }
    dims = dim(res)

  # check if inverse weights are given
  # this way if you pass Winv=NULL it will error out still
  if(is.null(Winv)){
    Winv = FALSE
    if(is.null(W)) W = rep(1,n)
  } else {
    W = Winv
    Winv = TRUE
  }

  # load mask
  if(class(mask)[1] !='niftiImage'){
    maskimg=as.character(mask)
    mask = RNifti::readNifti(maskimg)
  }

  # check that first input image and mask dimensions are the same
  ndims = length(dim(mask))
  if(any(dims[1:ndims] != dim(mask))  ){
    stop('images and mask dimensions do not match.\n')
  }

  # check that template and mask dimensions are the same
  if( !is.null(template)){
    if(class(template)[1]!='niftiImage'){
      temp = readNifti(as.character(template))
    } else {
      temp = template
    }
      dims = dim(temp)
      rm(temp)
      if(any(dims[1:ndims] != dim(mask))  ){
        stop('template image and mask dimensions (or pixel dimensions) do not match.\n')
      }
  }

  # load images
  res = t(apply(res, 4, function(x) x[mask!=0]))

  peind = which(!colnames(X) %in% colnames(Xred))
  df = length(peind)
  rdf = n - ncol(X)
 # if(df>1 & robust)
 #   stop('Robust covariance is only available for testing a single parameter.')

  if(is.character(W)){
    message('Weights are voxel-wise.\n')
    voxwts = TRUE
    W = do.call(abind::abind, list(RNifti::readNifti(W), along=4))
    W = t(apply(W, 4, function(x) x[mask!=0]))
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
    message('Running voxel-wise weighted linear models.\n')

    if(!robust){
      if(df==1){
        # cannot be parallelized due to memory use
        qrs = apply(W, 2, function(Wcol) qr(X * Wcol))
        # This should work in parallel
        if(.Platform$OS.type!='windows'){
          seX1 = sqrt(do.call(c, parallel::mclapply(qrs, function(qrval) chol2inv(qr.R(qrval))[peind,peind], mc.cores=mc.cores)))
        } else {
          seX1 = sqrt(do.call(c, lapply(qrs, function(qrval) chol2inv(qr.R(qrval))[peind,peind])))
        }
        # cannot be parallelized due to memory use. Also reshaping to a list to use mclapply takes longer.
        coef = num = do.call(cbind, lapply(1:ncol(res), function(ind) qr.coef(qrs[[ind]], res[,ind])[peind]) )
        # cannot be parallelized due to memory use
        res = do.call(rbind, lapply(1:ncol(res), function(ind) qr.resid(qrs[[ind]], res[,ind])) )
        rm(qrs, W)
      } else {
        # cannot be parallelized due to memory use. Also reshaping to a list to use mclapply takes longer.
        coef = num = do.call(cbind, lapply(1:ncol(res), function(ind) qr.coef(qrs[[ind]], res[,ind])[peind]) )
        num = do.call(c, lapply(1:ncol(res), function(ind) sum(qr.resid(qr(Xred * W[,ind]), res[,ind])^2)) )
        res = do.call(rbind, lapply(1:ncol(res), function(ind) qr.resid(qr(X * W[,ind]), res[,ind])) )
      }
    }

    if(robust){
      if(df==1){
      res = lapply(1:ncol(res), function(ind) lm(res[,ind] ~ -1 + I(X * W[,ind]), model=FALSE) )

      # get parameter estimates
      if(.Platform$OS.type!='windows'){
        coef = stat = do.call(cbind, parallel::mclapply(res, coefficients, mc.cores=mc.cores ))[peind,, drop=FALSE]
        message('Getting voxel-wise hat values.\n')
        h = do.call(rbind, parallel::mclapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}, mc.cores=mc.cores ))
        message('Getting voxel-wise residuals for covariate and outcome vectors.\n')
        res = do.call(rbind, parallel::mclapply(res, residuals, mc.cores=mc.cores))
      } else {
        coef = stat = do.call(rbind, lapply(res, coefficients ))[peind,,drop=FALSE]
        message('Getting voxel-wise hat values.\n')
        h = do.call(rbind, lapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}))
        message('Getting voxel-wise residuals for covariate and outcome vectors.\n')
        res = do.call(rbind, lapply(res, residuals))
      }

      if(!is.null(Xred)){
        X1res = do.call(rbind, lapply(1:nrow(res), function(ind) qr.resid(qr(Xred * W[,ind]), X1 * W[,ind])) )
      } else {
        # X1 is the intercept, Xred doesn't exist
        X1res = t(W) * X1
      }
      res = res * X1res /(1-h)
      A = rowSums(X1res^2)
      rm(h, X1res)
    } else {
      # compute qr decompositions
      qrs = lapply(1:ncol(res), function(ind) qr(X * W[,ind]) )
      # compute coefficients
      coef = do.call(cbind, lapply(1:ncol(res), function(ind) qr.coef(qrs[[ind]], res[,ind])[peind] ))
      # compute Q(v)
      temp = do.call(cbind, lapply(1:ncol(res), function(ind){r=qr.resid(qrs[[ind]], res[,ind]);
          h=rowSums(qr.Q(qrs[[ind]])^2); h = ifelse(h>=1, 1-eps, h)
          Q = W[,ind] * r/(1-h); Q }) )
      rm(qrs)
      # compute X_1^T W P^{X_0}. m1 X n X V array.
      # Depends on v here, but does not when weights are the same for all voxels
      system.time(
      test <- do.call(abind::abind, c(lapply(1:ncol(W), function(ind){qr.X0 = qr(Xred * W[,ind]); qr.resid(qr.X0, X1 * W[,ind]) }), along=3) ) )
      # now compute the 3d array that we need
      A = apply(test, 3, crossprod)
      # need this to simulate joint distribution
      res = sweep(test, c(1,3), temp, FUN = "*", check.margin=TRUE)
      rm(temp)
      # Compute Omega
      Omega = apply(res, 3, crossprod)
      # invert Omega
      Omega = apply(Omega, 2, function(x) chol2inv(chol(matrix(x, nrow=df, ncol=df))) )
      bA = do.call(cbind, lapply(1:ncol(A), function(ind) matrix(A[,ind], nrow=df, ncol=df) %*% coef[,1] ))
      stat = lapply(1:ncol(bA), function(ind) t(bA[,ind]) %*% matrix(Omega[,ind], nrow=df, ncol=df) %*% bA[,ind] )
      rm(bA, test)
    }
    }


  # else weights are the same for all voxels
  } else {
    QR = qr(X * W)
    coef = num = qr.coef(QR, res)[peind,,drop=FALSE]
    if(!robust){
      if(df==1){
        num = coef
        seX1 = sqrt(chol2inv(qr.R(QR))[peind,peind])
        res = t(qr.resid(QR, res))
	      rm(QR)
      } else {
        num = colSums(qr.resid(qr(Xred * W), res)^2)
        res = t(qr.resid(qr(X * W), res))
      }
    }

    if(robust){
      if(df==1){
      # message('Performing voxel regression.')
      # res = lm(res ~ -1 + I(X * W), model=FALSE)
      #
      # # get parameter estimates
      # coef = stat = coefficients(res)[peind,,drop=FALSE]
      #
      # # compute hat values
      # message('Computing hat values.')
      # h = rowSums(qr.Q(res$qr)^2)
      # h = ifelse(h>=1, 1-eps, h)
      #
      # # get residuals
      # message('Getting residuals.')
      # res = residuals(res)


      # qr approach
      message('Performing voxel regression.')
      qrX = qr(X * W)
      # PEs
      coef = stat = qr.coef(qrX, res)
      message('Computing hat values.')
      h = rowSums(qr.Q(qrX)^2)
      h = ifelse(h>=1, 1-eps, h)
      # residuals
      message('Getting residuals.')
      res = qr.resid(qrX, res)

      # residualize variable of interest to covariates
      # null statement is for if X is the intercept (Xred is null)
      X1res = if(!is.null(Xred)) qr.resid(qr(Xred * W), X1 * W) else X1 # Formula XX in the paper
      A = sum(X1res^2)
      # compute half of covariance of parameter of interest
      # divides by 1-h to use the HC3 version discussed by Long and Ervin
      # https://pdfs.semanticscholar.org/1526/72b624b44b12250363eee602554fe49ca782.pdf
      res = t(res *  (X1res/(1-h)))
      } else {
        message('blah')
      }
    }
  } # end if-else voxwts

  # compute statistical image
  if(!robust){
    message('Computing stat image.\n')
    stat = rowSums(res^2)
    # assume t-statistics if df==1
    if(df==1){
      stat = sqrt(stat/rdf)
      stat = num/stat / seX1
      # convert to z-statistics
      stat = qnorm(pt(stat, df=rdf))
    } else {
      num = num - stat
      stat = num/stat * rdf
      # convert to chisquared
      if(transform){
        stat = stat/df
        stat = qchisq(pf(stat, df1=df, df2=rdf), df=df)
      }
   }
  }

  if(robust){
    message('Computing robust stat image.\n')
    stat = stat*A/sqrt(rowSums(res^2))
    # Use T-to-Z transform
    if(transform) stat = qnorm(pt(stat, df=rdf))
  }

  if(!sqrtSigma) res=NULL


  # used later to indicated t-statistic
  if(df==1)
    df=0
  out = list(stat=stat, coef=coef, sqrtSigma=res, mask=mask, template=template, formulas=list(form, formred), robust=robust, df=df, rdf=rdf)
  class(out) = c('statMap', 'list')

  # if outdir is specified the stat and sqrtSigma images are saved in outdir
  # and mask tries to get saved as a character.
  if(!is.null(outdir)){
    files = write.statMap(out, outdir)
    out$stat = files$stat
    out$coef = files$coef
    out$sqrtSigma = files$sqrtSigma
    # if mask was a character then pass that forward instead if the niftiImage
    if(exists('maskimg'))
      out$mask = maskimg

  }
  return(out)
}
