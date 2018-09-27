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
#' @param robust Compute robust standard error estimates? Defaults to TRUE.
#'   Uses HC3 SE estimates from REF.
#' @param outdir If specified, output is saved as NIfTIs and statMap object is
#' saved as strings. This approach conserves memory, but has longer IO time.
#' Currently, not supported.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @importFrom abind abind
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical nifti object. If ncol(X) = ncol(Xred)+1, then this is a Z-statistic map, otherwise it is a chi^2-statistic map.}
#'   \item{res}{The 4d covariance object. This is a V by n matrix R, such that R \%*\% t(R) = Sigma.}
#' }
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @export
computeStats = function(images, form, formred, mask, data=NULL, W=rep(1, nrow(X)), Winv=NULL, template=NULL, formImages=NULL, robust=TRUE, outdir=NULL, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001

  # check if inverse weights are given
  if(!is.null(Winv)){
    W = Winv
    Winv = TRUE
  } else {
    Winv = FALSE
  }

  if(!is.matrix(form) & !is.matrix(formred)){
    X = getDesign(form, data)
    Xred = if(!is.null(formred)) getDesign(formred, data) else NULL
  } else {
    X = form
    Xred = formred
    form <- formred <- NULL
  }

  if(is.character(images)){
    n=length(images)
    if(nrow(X)!=n)
      stop('length(images) and nrow(X) must be the same.')
    res = do.call(abind::abind, list(RNifti::readNifti(images), along=4))
    } else {
      n = nrow(X)
      res = images
      rm(images)
    }

  # load mask
  if(is.character(mask)){
    maskimg=mask
    mask = RNifti::readNifti(mask)
  }

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
        if(.Platform$OS.type!='windows'){
          seX1 = sqrt(do.call(c, parallel::mclapply(qrs, function(qrval) chol2inv(qr.R(qrval))[peind,peind], mc.cores=mc.cores)))
        } else {
          seX1 = sqrt(do.call(c, lapply(qrs, function(qrval) chol2inv(qr.R(qrval))[peind,peind])))
        }
        # cannot be parallelized due to memory use. Also reshaping to a list to use mclapply takes longer.
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
      if(.Platform$OS.type!='windows'){
        stat = do.call(rbind, parallel::mclapply(res, coefficients, mc.cores=mc.cores ))[,peind]
        cat('Getting voxel-wise hat values.\n')
        h = do.call(rbind, parallel::mclapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}, mc.cores=mc.cores ))
        cat('Getting voxel-wise residuals for covariate and outcome vectors.\n')
        res = do.call(rbind, parallel::mclapply(res, residuals, mc.cores=mc.cores))
      } else {
        stat = do.call(rbind, lapply(res, coefficients ))[,peind]
        cat('Getting voxel-wise hat values.\n')
        h = do.call(rbind, lapply(res, function(r){ h=rowSums(qr.Q(r$qr)^2); h = ifelse(h>=1, 1-eps, h); h}))
        cat('Getting voxel-wise residuals for covariate and outcome vectors.\n')
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
      cat('Performing voxel regression.\n')
      res = lm(res ~ -1 + I(X * W), model=FALSE)

      # get parameter estimates
      stat=coefficients(res)[peind,]

      # compute hat values
      cat('Computing hat values.\n')
      h = rowSums(qr.Q(res$qr)^2)
      h = ifelse(h>=1, 1-eps, h)

      # get residuals
      cat('Getting residuals\n')
      res = residuals(res)

      # residualize variable of interest to covariates
      # null statement is for if X is the intercept (Xred is null)
      X1res = if(!is.null(Xred)) qr.resid(qr(Xred * W), X1 * W) else X1 # Formula XX in the paper
      A = sum(X1res^2)
      # compute half of covariance of parameter of interest
      # divides by 1-h to use the HC3 version discussed by Long and Ervin
      # https://pdfs.semanticscholar.org/1526/72b624b44b12250363eee602554fe49ca782.pdf
      res = t(res *  (X1res/(1-h)))
    }
  } # end if-else voxwts

  # compute statistical image
  if(!robust){
    cat('Computing stat image.\n')
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
    cat('Computing robust stat image.\n')
    stattemp = stat*A/sqrt(rowSums(res^2))
    # get niftiImage from mask
    stat = mask
    stat[ stat==1] = stattemp
  }


  # returns if requested
  out = list(stat=stat, sqrtSigma=res, mask=mask, template=template, formulas=list(form, formred), robust=robust)
  class(out) = c('statMap', 'list')

  # if outdir is specified the stat and sqrtSigma images are saved in outdir
  # and mask tries to get saved as a character.
  if(!is.null(outdir)){
    files = write.statMap(out, outdir)
    out$stat = files$stat
    out$sqrtSigma = files$sqrtSigma
    # if mask was a character then pass that forward instead if the niftiImage
    if(exists('maskimg'))
      out$mask = maskimg

  }
  return(out)
}
