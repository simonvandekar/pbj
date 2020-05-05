#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Spatial Extent Inference
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param statistic A user specified function that takes a RNift image object and computes a particular statistic of interest.
#' @param nboot Number of bootstrap samples to use.
#' @param rboot Function for generating random variables. See examples.
#' @param tboot Logical if an (approximate) t-bootstrap should be used. Currently, defaults to FALSE.
#' @param ... arguments passed to statistic function.
#'
#' @return Returns a list of length 3. The first element is the observed statistic value and the second is a list of the boostrap values.
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
pbjInference = function(statMap, statistic = function(image) max(c(image)), nboot=5000, rboot=stats::rnorm, tboot=FALSE, debug=FALSE, ...){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')


  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  rawstat = stat.statMap(statMap)
  template = statMap$template
  df = statMap$df
  rdf = statMap$rdf
  robust = statMap$robust


  if(df==0){
    sgnstat = sign(rawstat)
    stat = rawstat^2
    df=1; zerodf=TRUE
  } else {
    zerodf=FALSE
    stat = rawstat
  }

  obsstat = statistic(stat)

  sqrtSigma <- if(is.character(statMap$sqrtSigma)) {
    apply(readNifti(statMap$sqrtSigma), ndims+1, function(x) x[mask!=0])
  } else {
    statMap$sqrtSigma
  }
  rm(statMap)

  dims = dim(sqrtSigma)
  n = dims[2]
  V=dims[1]
  ndim = length(dims)
  if(is.null(rdf)) rdf=n

  if( length(dims)==2 ){
    # standardize these across n to be norm 1
    ssqs = sqrt(rowSums(sqrtSigma^2))
    # now Vxn
    if( any(ssqs != 1) ) sqrtSigma = t(sweep(sqrtSigma, 1, ssqs, '/'))
  } else {
    # This array should be standardized already
    # now nxVxm_1
    sqrtSigma = aperm(sqrtSigma, perm=c(2,1,3) )
  }
  # rearrange shape for ease below.
  dims = dim(sqrtSigma)

  #sqrtSigma <- as.big.matrix(sqrtSigma)
  boots = list()
  if(debug) statmaps = rep(list(NA), nboot)

  # if(.Platform$OS.type=='windows')
  # {

  bootDims = dim(rboot(n))
  bootLen = pmax(1,length(bootDims))
  if(is.null(bootDims)){
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
      if(length(dims)==2){
        if(is.null(bootDims)){
          boot = matrix(rboot(n*df), n, df)
        } else {
          boot = rboot(n)
        }
        # if df>1 this corresponds to an independence assumption among off-diagonal voxels
        # need to expand array to compute T-statistic
         if(tboot){
           if(df==1){
             statimg = sweep(sqrtSigma, 1, boot, FUN="*")
             statimg = apply(statimg, 2, function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
           } else {
             # independence approximation
             statimg = sweep(simplify2array(rep(list(sqrtSigma), df)), arrDims, boot, FUN="*")
             # standardize each voxel and normalized statistic
             statimg = apply(statimg, c(2,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
           }
         } else {
           statimg = t(sqrtSigma) %*% boot
         }
      } else {
        # assumes this is generating an n, nxV, or nxVxm_1 array
        boot = rboot(n)
        if(tboot){
          # dimensions are nxVxm_1
          statimg = sweep(sqrtSigma, arrDims, boot, FUN="*")
          # standardize each voxel and normalized statistic
          statimg = apply(statimg, c(2,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
        } else {
          statimg = colSums(sweep(sqrtSigma, arrDims, boot, FUN="*"), dims=1 )
        }
      }
      statimg = rowSums((statimg)^2)
      tmp[ mask!=0] = statimg
      boots[[i]] = statistic(tmp)
      setTxtProgressBar(pb, round(i/nboot,2))
    }
  close(pb)
  rm(sqrtSigma) # Free large big memory matrix object

  # add the stat max
  out = list(obsStat=obsstat, boots=boots)
  return(out)
}
