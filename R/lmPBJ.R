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
#' @param id Vector to identify measurements from the same observation.
#' @param data R data frame containing variables in form. If form and formred
#' are matrices then this can be NULL.
#' @param W Numeric vector of weights for regression model. Can be used to deweight noisy
#'  observations. Same as what should be passed to lm.
#' @param Winv Inverse weights for regression model. Inverse of W.
#' @param template Template image used for visualization.
#' @param formImages n X p matrix of images where n is the number of subjects and
#'  each column corresponds to an imaging covariate. Currently, not supported.
#' @param robust Logical, compute robust standard error estimates?
#' @param HC3 Logical, Uses HC3 SE estimates from Long and Ervin 2000? Defaults to TRUE.
#' @param transform character indicating type of transformation to use. "none", "t", "f", or "edgeworth" are currently accepted. Edgeworth is slow.
#' @param outdir If specified, output is saved as NIfTIs and statMap object is
#' saved as strings. This approach conserves memory, but has longer IO time.
#' Currently, not supported.
#' @param mc.cores Argument passed to mclapply for parallel things.
#' @param zeros Exclude voxels that have zeros? Zeros may exist due to differences in masking and
#' coverage or they may represent locations where the data took the value zero.
#' @keywords parametric bootstrap, statistical parametric map, semiparametric bootstrap
#' @return Returns a list with the following values:
#' \describe{
#'   \item{stat}{The statistical values where mask!=0. It is a chi^2-statistic map.}
#'   \item{coef}{A 4d niftiImage giving the parameter estimates for covariates only in the full model.}
#'   \item{sqrtSigma}{The covariance object used to sample from the joint distribution of the statistical image.}
#'   \item{mask}{The input mask.}
#'   \item{template}{The background template used for visualization.}
#'   \item{formulas}{A list containing the full and reduced models.}
#' }
#' stat=stat, sqrtSigma=res, mask=mask, template=template, formulas=list(form, formred), robust=robust, df=df, rdf=rdf
#' @importFrom stats coefficients lm pf pt qnorm qchisq residuals
#' @importFrom RNifti writeNifti readNifti
#' @importFrom parallel mclapply
#' @importFrom PDQutils papx_edgeworth
#' @export
#' @example inst/examples/lmPBJ.R
lmPBJ = function(images, form, formred=~1, mask, id=NULL, data=NULL, W=NULL,
Winv=NULL, template=NULL, formImages=NULL, robust=TRUE, transform=c('t', 'none', 'f', 'edgeworth'), outdir=NULL, zeros=FALSE, HC3=TRUE, mc.cores = getOption("mc.cores", 2L)){
  # hard coded epsilon for rounding errors in computing hat values
  eps=0.001
  transform = tolower(transform[1])

  X = getDesign(form, formred, data=data)
  Xred = X[['Xred']]
  df = X[['df']]
  images = images[X[['nas']] ]
  X = X[['X']]
  if(nrow(X) < nrow(data)){
    message(nrow(data)-nrow(X), ' observations deleted due to missingness.')
  }


  n = nrow(X)
  if(class(images[[1]])[1] != 'niftiImage'){
    if(n!=length(images))
      stop('length(images) and nrow(X) must be the same.')
    images = as.character(images)
    # removes white space after images if there is any
    images = gsub(" +$", "", images)
    Y = simplify2array(RNifti::readNifti(images))
  } else {
    Y = simplify2array(images)
    images=NULL
  }
  dims = dim(Y)

  if(is.character(W)){
    stop('Image valued weights are not supported.')
  }
  # check if inverse weights are given
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

  # inverse weights were passed
  if(Winv) W[W!=0] = 1/W[W!=0]
  X1 = X[,peind]
  W = sqrt(W)
  # this is a pointwise matrix multiplication if W was passed as images
  Y = Y * W

  # fit model to all image data
  QR = qr(X * W)
  coef = qr.coef(QR, Y)[peind,,drop=FALSE]
  res=qr.resid(QR, Y);

  X1res = qr.resid(qr(Xred * W), X1 * W)

  # standardize residuals and Y
  sigmas = sqrt(colSums(res^2)/rdf)
  #res = sweep(res, 2, sigmas, FUN = '/')
  #Y = sweep(Y, 2, sigmas, FUN = '/')

  if(!robust){
    AsqrtInv = backsolve(r=qr.R(qr(X1res)), x=diag(df) )
    sqrtSigma = crossprod(AsqrtInv, matrix(X1res, nrow=df, ncol=n, byrow=TRUE))
    # used to compute chi-squared statistic
    normedCoef = sweep(sqrtSigma %*% Y, 2, sigmas, FUN='/') # sweep((AsqrtInv%*% coef), 2, sigmas, FUN='/') #
    # In this special case only the residuals vary across voxels, so sqrtSigma can be obtained from the residuals
    sqrtSigma = list(res=res, X1res=as.matrix(X1res), QR=QR, XW=X*W, n=n, df=df, rdf=rdf, robust=robust, HC3=HC3, transform=transform)
    rm(AsqrtInv, Y, res, sigmas, X1res)
  } else {
    # first part of normedCoef
    normedCoef = colSums(sweep(simplify2array(rep(list(Y), df)), MARGIN = c(1,3), STATS = X1res, FUN = '*'), dims=1)
    if(HC3){
      h=rowSums(qr.Q(QR)^2); h = ifelse(h>=1, 1-eps, h)
      X1resQ = sweep(simplify2array(rep(list(res/(1-h)), df)),  c(1,3), X1res, '*')
    } else {
      # returns nXVXm_1 array
      X1resQ = sweep(simplify2array(rep(list(res), df)),  c(1,3), X1res, '*')
    }
    if(!is.null(id)){
      id = factor(id)
      IDmat = model.matrix(~-1+id)
      id = as.integer(id)
      X1resQ = array(apply(X1resQ, 3, function(mat) crossprod(IDmat, mat)), dim=c(ncol(IDmat), V, df))
    }
    # apply across voxels. returns V X m_1^2 array
    BsqrtInv = matrix(apply(X1resQ, 2, function(x){ backsolve(r=qr.R(qr(x)), x=diag(df)) }), nrow=df^2, ncol=V)
    #assign('BsqrtInvlmPBJ', BsqrtInv, envir = .GlobalEnv)
    # second part of normedCoef
    normedCoef = matrix(simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), normedCoef[ind,])) ), nrow=df)
    #assign('normedCoeflmPBJ', normedCoef, envir = .GlobalEnv)
    # Things needed to resample the robust statistics
    sqrtSigma = list(res=res, X1res=as.matrix(X1res), QR=QR, XW=X*W, n=n, df=df, rdf=rdf, robust=robust, HC3=HC3, transform=transform, id=id)
    rm(BsqrtInv, Y, res, X1resQ, X1res)
  }

  # use transform to compute chi-squared statistic
  normedCoef = switch(transform,
                      none=normedCoef,
                      f=normedCoef,
                      t={ qnorm(pt(normedCoef, df=rdf, log.p = TRUE ), log.p=TRUE )},
                      edgeworth={message('Computing edgeworth transform.')
                        matrix(qnorm(vpapx_edgeworth(stat=normedCoef, mu3=colSums(sqrtSigma^3, dims=1), mu4=colSums(sqrtSigma^4, dims=1) ) ), nrow=df)
                      })
  stat = colSums(normedCoef^2)
  if(transform=='f'){
    stat = qchisq(pf(stat/df, df1=df, df2=rdf, log.p = TRUE ), df=df, log.p=TRUE )
  }


  # used later to indicated t-statistic
  out = list(stat=stat, coef=coef, normedCoef=normedCoef, sqrtSigma=sqrtSigma, mask=mask, template=template, images=images, formulas=list(full=form, reduced=formred), data = get_all_vars(form, data = data))
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
