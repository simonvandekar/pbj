#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param statistic A user specified function that takes a RNift image object and computes a particular statistic of interest.
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. See examples.
#' @param method character method to use for bootstrap procedure.
#' @param ... arguments passed to statistic function.
#'
#' @return Returns a list of length 2. The first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @export
pbjInference = function(statMap, statistic = function(image) max(c(image)), nboot=5000, rboot=stats::rnorm, method=c('robust', 't', 'regular'), ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')


  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  ndims = length(dim(mask))
  rawstat = stat.statMap(statMap)
  template = statMap$template
  df = statMap$df
  rdf = statMap$rdf
  robust = statMap$robust
  stat = rawstat
  if(tolower(method[1]) == 'robust' & !robust) method = 't'

  obsstat = statistic(stat, ...)

  sqrtSigma <- if(is.character(statMap$sqrtSigma)) {
    apply(readNifti(statMap$sqrtSigma), ndims+1, function(x) x[mask!=0])
  } else {
    statMap$sqrtSigma
  }
  rm(statMap)

  dims = dim(sqrtSigma)
  n = dims[1]
  V=dims[2]
  ndim = length(dims)
  if(is.null(rdf)) rdf=n

  #sqrtSigma <- as.big.matrix(sqrtSigma)
  boots = list()

  # if(.Platform$OS.type=='windows')
  # {

  bootDims = dim(rboot(n))
  bootLen = pmax(1,length(bootDims))
  if(is.null(bootDims) & ndim==2){
    # the bootstrap is generating an n vector
    arrDims = c(1,3)
  } else {
    # assumes the bootstrap is generating an n X V array
    arrDims = 1:bootLen
  }

  if(!is.null(bootDims) & all(bootDims != dims[1:length(bootDims)] ) )
    stop('Dimension of bootstrap function does not match sqrtSigma')

    pb = txtProgressBar(style=3, title='Generating null distribution')
    tmp = mask
    for(i in 1:nboot){
      # bootstrap is univariate and sqrtSigma is 2D
      # assumes off-diagonal elements of spatial covariance are independent
      if(ndim==2){
        if(is.null(bootDims)){
          boot = matrix(rboot(n*df), n, df)
        } else {
          boot = rboot(n)
        }
        # if df>1 this corresponds to an independence assumption among off-diagonal voxels
        # need to expand array to compute T-statistic
        if(tolower(method[1])=='t'){
          if(df==1){
            statimg = sweep(sqrtSigma, 1, boot, FUN="*")
            statimg = apply(statimg, 2, function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
          } else {
            # off-diagonal spatially adjacent parameters are independent
            statimg = sweep(simplify2array(rep(list(sqrtSigma), df)), arrDims, boot, FUN="*")
            # standardize each voxel and normalized statistic
            statimg = apply(statimg, c(2,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
          }
        } else {
          statimg = crossprod(sqrtSigma, boot)
        }
      } else {
        # assumes this is generating an n or nxV array
        boot = rboot(n)
        # dimensions are nxVxm_1
        statimg = sweep(sqrtSigma, arrDims, boot, FUN="*")
        # only runs robust method if robust SE were used
        if(robust & tolower(method[1])=='robust'){
          BsqrtInv = matrix(apply(statimg, 2, function(x) pracma::sqrtm(crossprod(x))$Binv), nrow=df^2, ncol=V)
          statimg = t(simplify2array( lapply(1:V, function(ind) rowSums(tcrossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), matrix(statimg[,ind,], nrow=n, ncol=df))) ) ))
        } else if(tolower(method[1])=='regular'){
          statimg = colSums(statimg, dims=1 )
        } else {
          # standardize each voxel and normalized statistic
          statimg = apply(statimg, c(2,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
        }
      }
      statimg = rowSums((statimg)^2)
      tmp[ mask!=0] = statimg
      boots[[i]] = statistic(tmp, ...)
      setTxtProgressBar(pb, round(i/nboot,2))
    }
  close(pb)
  rm(sqrtSigma) # Free large big memory matrix object

  # add the stat max
  out = list(obsStat=obsstat, boots=boots)
  return(out)
}
