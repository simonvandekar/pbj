#' Performs (semi)Parametric Bootstrap Joint ((s)PBJ) Inference
#'
#' @param sqrtSigma list from statmap object
#' @param rboot a function that draws a bootstrapped sample. Should return an n vector. Defaults to Rademacher random variable.
#' @param null Is this a simulation under the null hypothesis?
#' @param method character, method to use for resampling procedure. Wild bootstrap, permutation, or nonparametric
#'
#' @return Returns vector of test statistics computed from the bootstrapped sample.
#' @export
#
pbjBoot = function(sqrtSigma, rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}, null=TRUE, method=c('wild', 'permutation', 'nonparametric')){
  method = tolower(method[1])
  eps=0.001
  V = ncol(sqrtSigma$res)
  id = sqrtSigma$id
  n = sqrtSigma$n
  df = sqrtSigma$df
  rdf = sqrtSigma$rdf
  HC3 = sqrtSigma$HC3
  robust = sqrtSigma$robust
  transform = sqrtSigma$transform
  if(HC3){
    h=rowSums(qr.Q(sqrtSigma$QR)^2); h = ifelse(h>=1, 1-eps, h)
    #h=rowSums(qr.Q(qr(sqrtSigma$X))^2); h = ifelse(h>=1, 1-eps, h)
  } else {
    h = rep(0, n)
  }


  if(robust){
    if(method=='wild'){#is.list(sqrtSigma)){ sqrtSigma should be a list here
      sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
    } else if (method=='permutation'){
      sqrtSigma$res = sqrtSigma$res[sample(n),]
    } else if (method=='nonparametric'){
      samp = sample(n, replace=TRUE)
      sqrtSigma$res = sweep(sqrtSigma$res[samp,], 1, sqrt(1-h[samp]), '/')
      sqrtSigma$X1res = sqrtSigma$X1res[samp,]
      sqrtSigma$XW = sqrtSigma$XW[samp,]
      sqrtSigma$QR = qr(sqrtSigma$XW)
    }
    # for bootstrapping under the alternative
    if(!null) sqrtSigma$res = sqrtSigma$XW %*% sqrtSigma$coef + sqrtSigma$res

    statimg = .Call("pbj_pbjBootRobustX", sqrtSigma$QR, sqrtSigma$res, sqrtSigma$X1res, id, h, df)
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
    # for bootstrapping under the alternative
    if(!null) sqrtSigma$res = sqrtSigma$XW %*% sqrtSigma$coef + sqrtSigma$res

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
                      f=statimg,
                      t={ qnorm(pt(statimg, df=rdf, log.p = TRUE ), log.p=TRUE )},
                      edgeworth={message('Computing edgeworth transform.')
                        matrix(qnorm(vpapx_edgeworth(stat=statimg, mu3=colSums(sqrtSigma$res^3, dims=1), mu4=colSums(sqrtSigma$res^4, dims=1) ) ), nrow=df)
                      })
  statimg = colSums(statimg^2)
  if(tolower(transform)=='f'){
    statimg = qchisq(pf(statimg/df, df1=df, df2=rdf, log.p = TRUE ), df=df, log.p=TRUE )
  }
  return(statimg)
}
