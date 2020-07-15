#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param boot a bootstrapped sample. Can be an n vector, V by n matrix, or n by df matrix.
#' @param V Number of voxels in analysis.
#' @param n Number of subjects.
#' @param df Numerator degrees of freedom.
#' @param method Method to use for resampling.
#' @param voxelwise logical indicating whether the data are voxelwise.
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, boot, V, n, df, method=c('robust', 't', 'conditional'), voxelwise=FALSE){
  method = tolower(method[1])
  if(!voxelwise){
    if(method=='robust'){#is.list(sqrtSigma)){ sqrtSigma should be a list here
      if( length(dim(boot))==0 ){ # dimension of bootstrap must be a vector of length n
        sqrtSigma$res = sweep(sqrtSigma$res, 1, boot, '*')
        BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(qr.resid(sqrtSigma$QR, sqrtSigma$res)), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
                                function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
        statimg = simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) ) ), higher=TRUE )
      } else {
        stop('Dimension of bootstrap sample is not correct for the method.')
      }
      } else if(method=='t'){
      if( all(dim(boot)==c(n,df)) ){ # dimension of bootstrap must be a matrix n X df of independent samples
        # in this case, off-diagonal spatially adjacent parameters are independent
        statimg = sweep(simplify2array(rep(list(sqrtSigma$res), df)), c(1,3), boot, FUN="*")
        # standardize each voxel and normalized statistic
        statimg = apply(statimg, c(1,3), function(x){ res = sum(x); res/sqrt(sum(x^2) -res^2/(length(x)-1) ) })
      } else {
        stop('Dimension of bootstrap sample is not correct for the method.')
      }
    } else if(method=='conditional'){
      statimg = crossprod(sqrtSigma$res, boot)
    }
  } else {
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
