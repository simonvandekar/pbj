#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param rboot a function that draws a bootstrapped sample. Can be an n vector, V by n matrix.
#' @param bootdim Dimension of rboot output
#' @param V Number of voxels in analysis.
#' @param n Number of subjects.
#' @param df Numerator degrees of freedom.
#' @param randomX Generate X in each bootstrap as though it is random as well.
#' @param robust Generate robust statistics?
#' @param method Method to use for resampling.
#' @param voxelwise logical indicating whether the data are voxelwise.
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, rboot, bootdim, V, n, df, randomX=FALSE, robust=TRUE, method=c('t', 'conditional', 'permutation'), voxelwise=FALSE){
  method = tolower(method[1])
  # !voxelwise
  if(!voxelwise){
    if(robust){
      if(method == 'conditional'){
        if( length(bootdim)==0 ){ # dimension of bootstrap must be a vector of length n
          boot = rboot(n)
        # could be made faster b/c BsqrtInv does not need to be computed in each bootstrap (though it is here)
        BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(qr.resid(sqrtSigma$QR, sqrtSigma$res)), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
                                function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
        sqrtSigma$res = sweep(sqrtSigma$res, 1, boot, '*')
        statimg = simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) ) ), higher=TRUE )
        } else {
          stop('Dimension of bootstrap sample is not correct for the method.')
        }
      } else{
        if(method=='t'){#is.list(sqrtSigma)){ sqrtSigma should be a list here
        if( length(bootdim)==0 ){ # dimension of bootstrap must be a vector of length n
          sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n), '*')
          if(randomX){
            sqrtSigma$X1res = sweep(sqrtSigma$X1res, 1, rboot(n), '*')
          }
        } else {
          stop('Dimension of bootstrap sample is not correct for the method.')
        }
      } else if (method=='permutation'){
          sqrtSigma$res = sqrtSigma$res[sample(n),]
      }
      # compute test statistic the regular way given the bootstrap/permuted sample
      BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(qr.resid(sqrtSigma$QR, sqrtSigma$res)), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
                              function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
      statimg = simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) ) ), higher=TRUE )
      }
    } else if(method=='t'){
      if( length(bootdim)==0 ){
        # dimension of bootstrap must be a matrix n X df of independent samples
        # in this case, off-diagonal spatially adjacent parameters are assumed to be independent
        boot = replicate(rboot(n), df)
      } else if(all(bootdim = c(n,df)) ){ # I don't think this will ever happen
        boot = rboot(n)
      } else {
        stop('Dimension of bootstrap sample is not correct for the method.')
      }
      statimg = sweep(simplify2array(rep(list(sqrtSigma$res), df)), c(1,3), boot, FUN="*")
      # standardize each voxel and normalized statistic
      statimg = apply(statimg, c(1,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
    } else if(method=='conditional'){
      boot = replicate(rboot(n), df)
      statimg = crossprod(boot, sqrtSigma$res)
    } else if(method=='permutation'){
      sqrtSigma$res = rboot(n)
      # compute t-statistic
    }

  } else { # voxelwise
    if(method=='robust'){
      stop('Voxelwise robust bootstrap not supported yet.')
    }
    if(method=='t'){
      stop('Voxelwise robust bootstrap not supported yet.')
    }
    if(method=='conditional'){
      stop('Voxelwise robust bootstrap not supported yet.')
    }
  }
  statimg = colSums(statimg^2)
}
