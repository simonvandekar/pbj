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
#' @param robust Logical, compute robust standard error estimates?
#' @param HC3 Logical, Uses HC3 SE estimates from Long and Ervin 2000? Defaults to TRUE.
#' @param transform character indicating type of transformation to use. "t" or "edgeworth." are currently accepted
#' (instead of niftiImage objects; defaults to FALSE).
#' @param outdir If specified, output is saved as NIfTIs and statMap object is
#' saved as strings. This approach conserves memory, but has longer IO time.
#' Currently, not supported.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @param zeros Exclude voxels that have zeros? Zeros may exist due to differences in masking and
#' coverage or they may represent locations where the data took the value zero.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical values where mask!=0. If ncol(X) = ncol(Xred)+1, then this is a Z-statistic map, otherwise it is a chi^2-statistic map.}
#'   \item{coef}{A 4d niftiImage giving the parameter estimates for covariates only in the full model.}
#'   \item{sqrtSigma}{The covariance object. This is a V by n by df matrix used for sampling from the statistical images joint distribution.}
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
#' @importFrom PDQutils papx_edgeworth
#' @export
lmPBJ = function(images, form, formred=~1, mask, data=NULL, W=NULL, Winv=NULL, template=NULL, formImages=NULL, robust=TRUE, transform=c('none', 't', 'edgeworth'), outdir=NULL, zeros=FALSE, HC3=TRUE, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001
  debug = getOption('pbj.debug', default=FALSE)

  X = getDesign(form, formred, data=data)
  Xred = X[['Xred']]
  df = X[['df']]
  X = X[['X']]

  if(class(images)[1] != 'niftiImage'){
    n=length(images)
    images = as.character(images)
    images = gsub(" +$", "", images)
    if(nrow(X)!=n)
      stop('length(images) and nrow(X) must be the same.')
    Y = simplify2array(RNifti::readNifti(images))
    } else {
      n = nrow(X)
      Y = images
      rm(images)
    }
    dims = dim(Y)

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
    stop('images and mask dimensions do not match.')
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
        stop('template image and mask dimensions (or pixel dimensions) do not match.')
      }
  }

  # load images
  if(zeros){
    # removes locations where there are any zeros
    mask = mask * c(apply(Y!=0, 1:ndims, all))
  }
  Y = t(apply(Y, (ndims+1), function(x) x[mask!=0]))
  V = ncol(Y)

  # assumes column names in X which aren't in Xred are of interest.
  peind = which(!colnames(X) %in% colnames(Xred))
  rdf = n - ncol(X) # this is true unless X is rank deficient

  if(is.character(W)){
    message('Weights are voxel-wise.')
    voxwts = TRUE
    W = simplify2array(RNifti::readNifti(W))
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
  Y = Y * W

  # fit model to all image data
  # if weights are voxel specific then design must also be treated separately
  if(voxwts){
    message('Running voxel-wise weighted linear models.')
    # cannot be parallelized due to memory use
    qrs = apply(W, 2, function(Wcol) qr(X * Wcol))# cannot be parallelized due to memory use. Also reshaping to a list to use mclapply takes longer.
    coef = simplify2array(  lapply(1:V, function(ind) qr.coef(qrs[[ind]], Y[,ind])[peind]) )

    if(!is.null(Xred)){
      X1res = simplify2array(lapply(1:V, function(ind) matrix(qr.resid(qr(Xred * W[,ind]), X1 * W[,ind]), nrow=n, ncol=df) ) )
    } else if(is.null(Xred) & df==1) {
      # X1 is the intercept, Xred doesn't exist nXm_1xV
      X1res = as.array(t(W) * X1, dim=c(nrow(W), 1, ncol(W) ))
    } else {
      stop('Degrees of freedom>1, but Xred is NULL.')
    }

    if(!robust){
      res = Y
      res = simplify2array( lapply(1:V, function(ind) qr.resid(qrs[[ind]], res[,ind])) )
      rm(qrs) # free memory
      # standardize residuals and Y
      sigmas = sqrt(colSums(res^2)/rdf)
      res = sweep(res, 2, sigmas, FUN = '/')
      Y = sweep(Y, 2, sigmas, FUN = '/')
      if(debug) assign('YlmPBJ', Y, envir = .GlobalEnv)
      AsqrtInv = matrix(apply(X1res, 3, function(x){ backsolve(r=qr.R(qr(x)), x=diag(df)) }),   nrow=df^2, ncol=V)
      sqrtSigma = simplify2array( lapply(1:V, function(ind) crossprod(matrix(AsqrtInv[,ind], nrow=df, ncol=df), matrix(X1res[,,ind], df, n, byrow=TRUE)) ) )
      # used to compute chi-squared statistic
      normedCoef = apply(sweep(sqrtSigma, MARGIN=c(2,3), Y, '*'), c(1,3), sum)
      # used for generating distribution of normedCoef
      sqrtSigma = sweep(sqrtSigma, MARGIN = c(2,3), res, '*' )
      sqrtSigma = aperm(sqrtSigma, c(2,3,1))
      rm(AsqrtInv, Y, res, sigmas, X1res)
    } else {
      # Only difference here is BsqrtInv instead of AsqrtInv and Q instead of res
      # compute Q(v)
      Q = simplify2array( lapply(1:V, function(ind){r=qr.resid(qrs[[ind]], Y[,ind]);
      if(HC3){
      h=rowSums(qr.Q(qrs[[ind]])^2); h = ifelse(h>=1, 1-eps, h)
      Q = r/(1-h)} else Q=r }) )
      rm(qrs) # free memory
      # first part of normedCoef
      normedCoef = colSums(sweep(X1res, MARGIN = c(1,3), Y, '*'), dims = 1)
      X1resQ = sweep(X1res, c(1,3), Q, '*')
      BsqrtInv = matrix(apply(X1resQ, 3, function(x){ backsolve(r=qr.R(qr(x)), x=diag(df)) }),   nrow=df^2, ncol=V)
      # second part of normedCoef
      normedCoef = matrix(simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), normedCoef[,ind])) ), nrow=df)
      sqrtSigma = X1resQ
      sqrtSigma = aperm(sqrtSigma, c(1,3,2))
      rm(BsqrtInv, Y, Q, X1resQ, X1res)
    }
    # weights are not voxel specific
  } else {

    message('Running weighted linear models.')
    QR = qr(X * W)
    coef = qr.coef(QR, Y)[peind,,drop=FALSE]
    res=qr.resid(QR, Y);

    if(!is.null(Xred)){
      X1res = qr.resid(qr(Xred * W), X1 * W)
    } else if(is.null(Xred) & df==1) {
      # X1 is the intercept, Xred doesn't exist nXm_1xV
      X1res = as.matrix(X1 * W)
    } else {
      stop('Degrees of freedom>1, but Xred is NULL.')
    }

    if(!robust){
      # standardize residuals and Y
      sigmas = sqrt(colSums(res^2)/rdf)
      res = sweep(res, 2, sigmas, FUN = '/')
      Y = sweep(Y, 2, sigmas, FUN = '/')
      if(debug) assign('YlmPBJ', Y, envir = .GlobalEnv)
      AsqrtInv = backsolve(r=qr.R(qr(X1res)), x=diag(df) )
      sqrtSigma = crossprod(AsqrtInv, matrix(X1res, nrow=df, ncol=n, byrow=TRUE))
      # used to compute chi-squared statistic
      normedCoef = sqrtSigma %*% Y # sweep((AsqrtInv%*% coef), 2, sigmas, FUN='/') #
      # In this special case only the residuals vary across voxels, so sqrtSigma can be obtained from the residuals
      sqrtSigma = list(res=res, X1res=X1res, QR=QR, XW=X*W, n=n, df=df, rdf=rdf, Y=Y)
      rm(AsqrtInv, Y, res, sigmas, X1res)
    } else {
      # standardize residuals and Y
      sigmas = sqrt(colSums(res^2)/rdf)
      res = sweep(res, 2, sigmas, FUN = '/')
      Y = sweep(Y, 2, sigmas, FUN = '/')
      if(debug) assign('YlmPBJ', Y, envir = .GlobalEnv)
      # first part of normedCoef
      normedCoef = colSums(sweep(simplify2array(rep(list(Y), df)), MARGIN = c(1,3), STATS = X1res, FUN = '*'), dims=1)
      if(HC3){
        h=rowSums(qr.Q(QR)^2); h = ifelse(h>=1, 1-eps, h)
        X1resQ = sweep(simplify2array(rep(list(res/(1-h)), df)),  c(1,3), X1res, '*')
      } else {
        # returns nXVXm_1 array
        X1resQ = sweep(simplify2array(rep(list(res), df)),  c(1,3), X1res, '*')
      }
      # apply across voxels. returns V X m_1^2 array
      BsqrtInv = matrix(apply(X1resQ, 2, function(x){ backsolve(r=qr.R(qr(x)), x=diag(df)) }), nrow=df^2, ncol=V)
      #assign('BsqrtInvlmPBJ', BsqrtInv, envir = .GlobalEnv)
      # second part of normedCoef
      normedCoef = matrix(simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), normedCoef[ind,])) ), nrow=df)
      #assign('normedCoeflmPBJ', normedCoef, envir = .GlobalEnv)
      # Things needed to resample the robust statistics
      sqrtSigma = list(res=res, X1res=X1res, QR=QR, XW=X*W, n=n, df=df, rdf=rdf)
      rm(BsqrtInv, Y, res, X1resQ, X1res)
    }
  }

  # use transform to compute chi-squared statistic
  normedCoef = switch(tolower(transform[1]),
         none=normedCoef,
         t={ qnorm(pt(normedCoef, df=rdf ) )},
         edgeworth={message('Computing edgeworth transform.')
           matrix(qnorm(vpapx_edgeworth(stat=normedCoef, mu3=colSums(sqrtSigma^3, dims=1), mu4=colSums(sqrtSigma^4, dims=1) ) ), nrow=df)
         })
  stat = colSums(normedCoef^2)


  # used later to indicated t-statistic
  out = list(stat=stat, coef=coef, normedCoef=normedCoef, sqrtSigma=sqrtSigma, mask=mask, template=template, formulas=list(form, formred), robust=robust, df=df, rdf=rdf, HC3=HC3, transform=transform)
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
