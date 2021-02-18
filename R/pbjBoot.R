#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param rboot a function that draws a bootstrapped sample. Should return an n vector. Defaults to Rademacher random variable.
#' @param bootdim Dimension of rboot output
#' @param robust Generate robust statistics?
#' @param transform Apply a quantile transformation to the test statistics to improve normal approximation.
#' @param method Method to use for resampling.
#' @param voxelwise logical indicating whether the data are voxelwise.
#' @param HC3 logical, was HC3 estimator used?
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, bootdim, method=c('t', 'independence', 'permutation', 'conditional', 'nonparametric'), voxelwise=FALSE, HC3=TRUE, robust=TRUE, transform=c('none', 't')){
  method = tolower(method[1])
  eps=0.001
  V = ncol(sqrtSigma$res)
  n = sqrtSigma$n
  df = sqrtSigma$df
  rdf = sqrtSigma$rdf
  # !voxelwise
  if(!voxelwise){
    if(HC3){
      h=rowSums(qr.Q(sqrtSigma$QR)^2); h = ifelse(h>=1, 1-eps, h)
      #h=rowSums(qr.Q(qr(sqrtSigma$X))^2); h = ifelse(h>=1, 1-eps, h)
    } else {
      h = rep(0, n)
    }
    if(robust){
      if(method == 'conditional'){
        boot = rboot(n)
        # could be made faster b/c BsqrtInv does not need to be computed in each bootstrap (though it is here)
        BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(sweep(qr.resid(sqrtSigma$QR, sqrtSigma$res), 1, 1-h, '/')), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
                                function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
        sqrtSigma$res = sweep(sqrtSigma$res, 1, boot/sqrt(1-h), '*')
        statimg = matrix(simplify2array( lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) ) ), higher=TRUE ), nrow=df, ncol=V)
      } else if(method=='independence'){
        # n x df array
        boot = replicate(df, rboot(n))
        # Divide by 1-h for debiased residuals. result is n x V x df array
        statimg = sweep(simplify2array(rep(list(sweep(sqrtSigma$res, 1, 1-h, '/')), df), higher = TRUE), c(1,3), sqrtSigma$X1res * boot, '*')
        # scale to unit variance and compute standardized parameter estimates.
        statimg = matrix(colMeans(apply(statimg, c(2,3), scale), dims = 2), nrow=df, ncol=V)
      } else{
        if(method=='t'){#is.list(sqrtSigma)){ sqrtSigma should be a list here
          sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
          #sigmas = sqrt(colSums(qr.resid(sqrtSigma$QR, sqrtSigma$res)^2)/rdf)
          #sqrtSigma$res = sweep(sqrtSigma$res, 2, sigmas, FUN = '/')
        } else if (method=='permutation'){
          sqrtSigma$res = sqrtSigma$res[sample(n),]
          #sqrtSigma$res = sqrtSigma$res[1:n,]
        } else if (method=='nonparametric'){
          samp = sample(n, replace=TRUE)
          sqrtSigma$res = sweep(sqrtSigma$res[samp,], 1, sqrt(1-h[samp]), '/')
          sqrtSigma$X1res = sqrtSigma$X1res[samp,]
          sqrtSigma$XW = sqrtSigma$XW[samp,]
          sqrtSigma$QR = qr(sqrtSigma$XW)
        }
        #else if(method=='robustpermutation'){
        #sqrtSigma$res = sqrt(abs(sqrtSigma$res[sample(n),])) * sign(sqrtSigma$res) * sqrt(abs(sqrtSigma$res))
        #sqrtSigma$res = sqrtSigma$res[sample(n),] * abs(sqrtSigma$res)
        #}
        # compute test statistic the regular way given the bootstrap/permuted sample
        #BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(qr.resid(sqrtSigma$QR, sqrtSigma$res)), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
        #                        function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
        BsqrtInv = matrix(apply(sweep(simplify2array(rep(list(sweep(qr.resid(sqrtSigma$QR, sqrtSigma$res), 1, 1-h, '/')), df)), c(1,3), sqrtSigma$X1res, '*'), 2,
                                function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
        statimg = matrix(simplify2array(lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) ) ), higher=TRUE ), nrow=df, ncol=V)
      }
    } else {
      if(method=='conditional'){
        # use this method if transform = 't'
        sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
      } else {
        if(method=='t'){
          sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
        } else if(method=='permutation'){
          sqrtSigma$res = sqrtSigma$res[sample(n), ]
        } else if (method=='nonparametric'){
          samp = sample(n, replace=TRUE)
          sqrtSigma$res = sweep(sqrtSigma$res[samp,], 1, sqrt(1-h[samp]), '/')
          sqrtSigma$X1res = sqrtSigma$X1res[samp,]
          sqrtSigma$XW = sqrtSigma$XW[samp,]
          sqrtSigma$QR = qr(sqrtSigma$XW)
        }
        sigmas = sqrt(colSums(qr.resid(sqrtSigma$QR, sqrtSigma$res)^2)/(rdf))
        sqrtSigma$res = sweep(sqrtSigma$res, 2, sigmas, FUN = '/')
      }
      AsqrtInv = backsolve(r=qr.R(qr(sqrtSigma$X1res)), x=diag(df) )
      statimg = crossprod(AsqrtInv, matrix(sqrtSigma$X1res, nrow=df, ncol=n, byrow=TRUE))
      # used to compute chi-squared statistic
      statimg = statimg %*% sqrtSigma$res
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
  statimg = switch(tolower(transform[1]),
                   none=statimg,
                   t={ qnorm(pt(statimg, df=rdf ) )})
  statimg = colSums(statimg^2)
}
