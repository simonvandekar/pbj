test_that("pbj_pbjBootRobustX performs correctly when: method = 'wild', robust = TRUE, HC3 = TRUE, transform = 'none'", {
  # setting up
  pain <- pain21::pain21()
  pain$data$group <- factor(sample(1:4, size = nrow(pain$data), replace=TRUE))
  pain$data$x <- rnorm(nrow(pain$data))
  pain$data$Winv <- runif(nrow(pain$data))
  imgs <- simplify2array(RNifti::readNifti(pain$data$images))
  Winvs <- simplify2array(RNifti::readNifti(pain$data$varimages))
  mask <- RNifti::readNifti(pain$mask) * c(apply(imgs!=0, 1:3, all))

  # get one voxel
  testvox <- which(mask==1, arr.ind = TRUE)[1:2,]
  mask[,,] <- 0
  mask[testvox] <- 1

  pain$data$y <- imgs[testvox[1,1], testvox[1,2], testvox[1,3], ]
  pain$data$Winv.img <- Winvs[testvox[1,1], testvox[1,2], testvox[1,3], ]

  outdir <- tempfile("dir")
  dir.create(outdir)
  statMap <- lmPBJ(pain$data$images, form = ~ splines::ns(x, df = 3), formred = ~ 1, mask = mask,
                   template = pain$template, data = pain$data,
                   Winv = pain$data$Winv, zeros = TRUE, transform = 'none',
                   robust = TRUE, HC3 = TRUE, outdir = outdir)

  sqrtSigma <- readRDS(statMap$sqrtSigma)
  mask <- statMap$mask
  dims <- dim(sqrtSigma$res)
  rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}

  # R implementation
  set.seed(1234)
  eps <- 0.001
  V <- ncol(sqrtSigma$res)
  n <- sqrtSigma$n
  df <- sqrtSigma$df
  h <- rowSums(qr.Q(sqrtSigma$QR)^2)
  h <- ifelse(h >= 1, 1 - eps, h)

  res <- sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
  BsqrtInv <- matrix(apply(sweep(simplify2array(rep(list(sweep(qr.resid(sqrtSigma$QR, res), 1, 1-h, '/')), df)), c(1,3), sqrtSigma$X1res, '*'), 2, function(x){ backsolve(r=qr.R(qr(x)), x=diag(ncol(x))) }), nrow=df^2, ncol=V)
  statimg <- matrix(simplify2array(lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, res[,ind]) ) ), higher=TRUE ), nrow=df, ncol=V)
  expected <- colSums(statimg^2)

  # C implementation
  set.seed(1234)
  actual <- pbjBoot(sqrtSigma, method = 'wild')

  expect(all.equal(actual, expected), "C implementation didn't match")
  expect(identical(sqrtSigma, readRDS(statMap$sqrtSigma)), "sqrtSigma was modified in-place")
})



test_that("pbj_pbjBootRobustX performs correctly when: method = 'wild', robust = TRUE, HC3 = TRUE, id!=NULL, transform = 'none'", {
  # setting up
  pain <- pain21::pain21()
  pain$data$group <- factor(sample(1:4, size = nrow(pain$data), replace=TRUE))
  pain$data$x <- rnorm(nrow(pain$data))
  pain$data$Winv <- runif(nrow(pain$data))
  # creates a fake ID variable
  pain$data$ID = c(rep(1:10, each=2), 11)
  imgs <- simplify2array(RNifti::readNifti(pain$data$images))
  Winvs <- simplify2array(RNifti::readNifti(pain$data$varimages))
  mask <- RNifti::readNifti(pain$mask) * c(apply(imgs!=0, 1:3, all))

  # get one voxel
  testvox <- which(mask==1, arr.ind = TRUE)[1:2,]
  mask[,,] <- 0
  mask[testvox] <- 1

  pain$data$y <- imgs[testvox[1,1], testvox[1,2], testvox[1,3], ]
  pain$data$Winv.img <- Winvs[testvox[1,1], testvox[1,2], testvox[1,3], ]

  outdir <- tempfile("dir")
  dir.create(outdir)
  statMap <- lmPBJ(pain$data$images, form = ~ 1, formred = NULL, mask = mask,
                   template = pain$template, data = pain$data, id = pain$data$ID,
                   Winv = pain$data$Winv, zeros = TRUE, transform = 'none',
                   robust = TRUE, HC3 = TRUE, outdir = outdir)

  sqrtSigma <- readRDS(statMap$sqrtSigma)
  mask <- statMap$mask
  dims <- dim(sqrtSigma$res)
  rboot=function(n){ (2*stats::rbinom(n, size=1, prob=0.5)-1)}

  # R implementation
  set.seed(1234)
  eps <- 0.001
  V <- ncol(sqrtSigma$res)
  n <- sqrtSigma$n
  df <- sqrtSigma$df
  h <- rowSums(qr.Q(sqrtSigma$QR)^2)
  h <- ifelse(h >= 1, 1 - eps, h)
  id = pain$data$ID

  res <- sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
  X1resQ = sweep(simplify2array(rep(list(sweep(qr.resid(sqrtSigma$QR, res), 1, 1-h, '/')), df)), c(1,3), sqrtSigma$X1res, '*')
  if(!is.null(id)){
    id = factor(id)
    IDmat = model.matrix(~-1+id)
    id = as.numeric(id)
    X1resQ = array(apply(X1resQ, 3, function(mat) crossprod(IDmat, mat)), dim=c(ncol(IDmat), V, df))
  }
  # apply across voxels. returns V X m_1^2 array
  BsqrtInv = matrix(apply(X1resQ, 2, function(x){ backsolve(r=qr.R(qr(x)), x=diag(df)) }), nrow=df^2, ncol=V)
  statimg <- matrix(simplify2array(lapply(1:V, function(ind) crossprod(matrix(BsqrtInv[,ind], nrow=df, ncol=df), crossprod(sqrtSigma$X1res, res[,ind]) ) ), higher=TRUE ), nrow=df, ncol=V)
  expected <- colSums(statimg^2)

  #X1resQ = matrix(as.numeric(strsplit("0.141684 -5.573603 -0.648924 1.048352 6.923937 -0.360543 -0.261023 -2.631793 0.002282 0.264740 0.253521 -0.102001 -6.513960 -0.252866 0.557416 6.684720 -0.393991 -0.449779 -2.031558 0.966052 0.280116 0.276790", split=" ")[[1]]), nrow = 11)

  # C implementation
  set.seed(1234)
  actual <- pbjBoot(sqrtSigma, method = 'wild')

  expect(all.equal(actual, expected), "C implementation didn't match")
  expect(identical(sqrtSigma, readRDS(statMap$sqrtSigma)), "sqrtSigma was modified in-place")
})
