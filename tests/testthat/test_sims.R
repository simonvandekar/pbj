# sim test
# test that these commands work

library(pain21)
library(testthat)
library(lmtest)
library(sandwich)
library(splines)
#devtools::load_all('./')

# setting up
set.seed(1234)
pain = pain21::pain21()
pain$data$group = factor(sample(1:4, size = nrow(pain$data), replace=TRUE))
pain$data$x = rnorm(nrow(pain$data))
pain$data$Winv = runif(nrow(pain$data))
#debug(lmPBJ)
# test by comparing one voxel to results obtained by lmtest and sandwich packages
imgs = simplify2array(RNifti::readNifti(pain$data$images))
Winvs = simplify2array(RNifti::readNifti(pain$data$varimages))
mask = RNifti::readNifti(pain$mask) * c(apply(imgs!=0, 1:3, all))

# get two voxels
testvox = which(mask==1, arr.ind = TRUE)[1:2,]
mask[,,] = 0
mask[testvox] = 1

# Statistic function to get objects for pbjInference
simStats = function(image, mask, thrs){
  c(maximum = max(c(image)), pbj::cluster(image, mask, thrs))
}

# simfunc should contain a data argument, which is defined within runSim
# Other arguments are identical across simulation runs.
simFunc = function(lmfull, lmred, mask, data, nboot, cfts){
  # generate fake covariates
  data$fake_group = factor(ceiling(ppoints(nrow(data))*4 ) )
  data$fake_covariate1 = rnorm(nrow(data))

  # compute Statistical maps
  # Residuals from null model to use for permutation procedures
  nullResids = lmPBJ(data$images, form=lmred, formred=~ 1, mask=mask, data=data, transform = 'none', robust=FALSE )$sqrtSigma
  perm = function(n){
    nullResids[sample(1:n),]
  }
  signedResids = sign(nullResids)
  exchPerm = function(n){
    signedResids[sample(1:n),]
  }

  # t transform, robust, estimate covariance
  tRobustStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 't' )
  # t transform, classical, estimate covariance
  tStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 't', robust=FALSE)
  # Doesn't scale residuals by hat matrix diagonal
  tPermStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 't', HC3=FALSE )
  # no transform, robust, estimate covariance
  #robustStatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform ='none')
  # no transform, classical, estimate covariance
  plainstatmap = lmPBJ(data$images, form=lmfull, formred=lmred, mask=mask, data=data, transform = 'none', robust=FALSE)
  resids = plainstatmap$sqrtSigma

  statmaps = c('tRobustStatmap', 'tStatmap', 'tPermStatmap')
  out = list()
  # doesn't matter which statmap we use here
  thrs = (cfts^2*plainstatmap$rdf) + plainstatmap$df
  # Apply each of the sampling methods
  for(statmapname in statmaps){

    ### BOOTSTRAP METHODS
    statmap = get(statmapname)
    # normal bootstrap
    #pbjNorm = getBoots(pbjSEI(statmap, nboot = nboot, cfts.s = cfts))
    if(statmapname %in% c('tRobustStatmap', 'tStatmap')){
      pbjNormT = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, method='t')
      # Rademacher bootstrap
      pbjRadT = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, rboot = function(n){ 2*rbinom(n, size=1, prob=0.5)-1}, method='t')
    }
    if(statmapname=='tRobustStatmap'){
      pbjRadRobust = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, rboot = function(n){ 2*rbinom(n, size=1, prob=0.5)-1}, method='robust')
      pbjNormRobust = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, method='robust')
      # exchangeable permutation
      # scale out sign of Q from sqrtSigma -- sign comes from "permExch" function
      statmap$sqrtSigma = sweep(statmap$sqrtSigma, MARGIN = c(1,2), STATS = sign(resids), FUN="*")
      pbjPermExch = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, rboot = exchPerm, method='robust')
    }
    if(statmapname == 'tStatmap'){
      pbjNorm = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, method='regular')
      pbjRad = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, rboot = function(n){ 2*rbinom(n, size=1, prob=0.5)-1}, method='regular')
    }

    if(statmapname=='tPermStatmap'){
      # now scale out residuals from sqrtSigma -- scale comes from "perm" function
      statmap$sqrtSigma = sweep(statmap$sqrtSigma, MARGIN = c(1,2), STATS = resids, FUN="/")
      pbjPerm = pbjInference(statmap, nboot = nboot, thr = thrs, mask=statmap$mask, statistic=simStats, rboot = perm, method='robust')
    }
    # collect output
    PBJnames = grep('^pbj', ls(), value=TRUE)
    allnames = paste(statmapname, PBJnames, sep='_')
    out[allnames] = lapply(PBJnames, get, pos = environment())
    rm(PBJnames)

    ### REPEAT ALL WITH INDEPENDENCE SPATIAL COVARIANCE ASSUMPTION
    # nonrobust methods won't be different, because covariance is same for all statistics.
  }
  return(out)
}

test = simFunc(lmfull = ~ fake_group + fake_covariate1, lmred = ~ fake_covariate1,
        mask=mask, data=pain$data, nboot=10, cfts=c(0.1, 0.2))

