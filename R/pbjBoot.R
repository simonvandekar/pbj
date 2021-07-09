#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param rboot a function that draws a bootstrapped sample. Should return an n vector. Defaults to Rademacher random variable.
#' @param bootdim Dimension of rboot output
#' @param robust Generate robust statistics?
#' @param transform Apply a quantile transformation to the test statistics to improve normal approximation.
#' @param method character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#' @param HC3 logical, was HC3 estimator used?
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, bootdim, method=c('wild', 'permutation', 'nonparametric'), HC3=TRUE, robust=TRUE, transform=c('none', 't')){
  method = tolower(method[1])
  eps=0.001
  V = ncol(sqrtSigma$res)
  n = sqrtSigma$n
  df = sqrtSigma$df
  rdf = sqrtSigma$rdf
    if(HC3){
      h=rowSums(qr.Q(sqrtSigma$QR)^2); h = ifelse(h>=1, 1-eps, h)
      #h=rowSums(qr.Q(qr(sqrtSigma$X))^2); h = ifelse(h>=1, 1-eps, h)
    } else {
      h = rep(0, n)
    }
    if(robust){
        if(method=='wild'){#is.list(sqrtSigma)){ sqrtSigma should be a list here
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

        statimg = .Call("pbj_pbjBootRobustX", sqrtSigma$QR, sqrtSigma$res, sqrtSigma$X1res, h, df)
    } else {
        if(method=='wild'){
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
      # this could be performed outside of the bootstrap function
      AsqrtInv = backsolve(r=qr.R(qr(sqrtSigma$X1res)), x=diag(df) )
      statimg = crossprod(AsqrtInv, matrix(sqrtSigma$X1res, nrow=df, ncol=n, byrow=TRUE))
      # used to compute chi-squared statistic
      statimg = statimg %*% sqrtSigma$res
    }


  statimg = switch(tolower(transform[1]),
                   none=statimg,
                   t={ qnorm(pt(statimg, df=rdf ) )})
  statimg = colSums(statimg^2)
}
