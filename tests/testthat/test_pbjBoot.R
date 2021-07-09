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
  statMap <- lmPBJ(pain$data$images, form = ~ 1, formred = NULL, mask = mask,
                   template = pain$template, data = pain$data,
                   Winv = pain$data$Winv, zeros = TRUE, transform = 'none',
                   robust = TRUE, HC3 = TRUE, outdir = outdir)

  rboot <- function(n) { (2*stats::rbinom(n, size=1, prob=0.5)-1) }
  sqrtSigma <- readRDS(statMap$sqrtSigma)
  mask <- statMap$mask
  dims <- dim(sqrtSigma$res)
  bootdim <- dim(rboot(nrow(sqrtSigma$res)))

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
  actual <- pbjBoot(sqrtSigma, rboot, bootdim, robust = TRUE, method = 'wild',
                    HC3 = TRUE, transform = 'none')

  expect(identical(actual, expected), "C implementation didn't match")
  expect(identical(sqrtSigma, readRDS(statMap$sqrtSigma)), "sqrtSigma was modified in-place")
})
