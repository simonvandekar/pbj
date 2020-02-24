# Tests that pbj output matches univariate methods

# first install and restart
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

# get one voxel
testvox = which(mask==1, arr.ind = TRUE)[1:2,]
mask[,,] = 0
mask[testvox] = 1

# use R functions from sandwich and lmtest
pain$data$y = imgs[testvox[1,1], testvox[1,2], testvox[1,3], ]
pain$data$Winv.img = Winvs[testvox[1,1], testvox[1,2], testvox[1,3], ]


# a series of tests to see if my code matches standard R output
tol = 10^(-6)

#### scalar weights df=2
test_that("Output from PBJ with df=2 and scalar weights matches output from lmtest and sandwich packages.", {
  model = lm(y ~ group, data=pain$data, weights =  1/pain$data$Winv)
  model.red = lm(y ~ 1, data=pain$data, weights=1/pain$data$Winv)
  waldtestres = lmtest::waldtest(model, model.red, test='F', vcov=sandwich::vcovHC)
  # pbj methods
  statmap <- lmPBJ(pain$data$images, form = ~ group,
                   formred = ~ 1, mask = mask,
                   template=pain$template, data = pain$data,
                   Winv = pain$data$Winv, zeros=TRUE)
  expect_equal(statmap$coef[,1], coefficients(model)[-1], tolerance=tol )
  expect_equal(qf(pchisq(statmap$stat, df=statmap$df), df1 = statmap$df, df2=statmap$rdf)[1], waldtestres$F[2], tolerance=tol)
})


#### voxel-wise weights df=2
test_that("Output from PBJ with df=2 and image weights matches output from lmtest and sandwich packages.", {
  model = lm(y ~ group, data=pain$data, weights =  1/pain$data$Winv.img)
  model.red = lm(y ~ 1, data=pain$data, weights=1/pain$data$Winv.img)
  waldtestres = lmtest::waldtest(model, model.red, test='F', vcov=sandwich::vcovHC)
  statmap <- lmPBJ(pain$data$images, form = ~ group,
                   formred = ~ 1, mask = mask,
                   template=pain$template, data = pain$data,
                   Winv = pain$data$varimages, zeros=TRUE)
  expect_equal(statmap$coef[,1], coefficients(model)[-1], tolerance=tol )
  expect_equal(qf(pchisq(statmap$stat, df=statmap$df), df1 = statmap$df, df2=statmap$rdf)[1], waldtestres$F[2], tolerance=tol)
})
# scalar weights df=1
test_that("Output from PBJ with df=1 and scalar weights matches output from lmtest and sandwich packages.", {
  model = lm(y ~ x, data=pain$data, weights =  1/pain$data$Winv)
  model.red = lm(y ~ 1, data=pain$data, weights=1/pain$data$Winv)
  waldtestres = lmtest::waldtest(model, model.red, test='F', vcov=sandwich::vcovHC)
  statmap <- lmPBJ(pain$data$images, form = ~ x,
                   formred = ~ 1, mask = mask,
                   template=pain$template, data = pain$data,
                   Winv = pain$data$Winv, zeros=TRUE)
  expect_equal(statmap$coef[,1], coefficients(model)[-1], tolerance=tol )
  expect_equal(qf(pchisq(statmap$stat^2, df=1), df1 = 1, df2=statmap$rdf)[1], waldtestres$F[2], tolerance=tol)
})

# voxel-wise weights df=1
test_that("Output from PBJ with df=1 and image weights matches output from lmtest and sandwich packages.", {
  model = lm(y ~ x, data=pain$data, weights =  1/pain$data$Winv.img)
  model.red = lm(y ~ 1, data=pain$data, weights=1/pain$data$Winv.img)
  waldtestres = lmtest::waldtest(model, model.red, test='F', vcov=sandwich::vcovHC)
  statmap <- lmPBJ(pain$data$images, form = ~ x,
                   formred = ~ 1, mask =mask,
                   template=pain$template, data = pain$data,
                   Winv = pain$data$varimages, zeros=TRUE)
  expect_equal(statmap$coef[,1], coefficients(model)[-1], tolerance=tol )
  expect_equal(qf(pchisq(statmap$stat^2, df=1), df1 = 1, df2=statmap$rdf)[1], waldtestres$F[2], tolerance=tol)
} )


test_that("Output from PBJ with nonlinear test and scalar weights matches output from lmtest and sandwich packages.", {
  models = getDesign(~ ns(x, df = 4), ~ x, data=pain$data)
  pain$data[, c('x', 'u1', 'u2', 'u3') ] = models$X[,-1]
  model = lm(y ~ x + u1 + u2 + u3, data=pain$data, weights =  1/pain$data$Winv)
  model.red = lm(y ~ x, data=pain$data, weights=1/pain$data$Winv)
  waldtestres = lmtest::waldtest(model, model.red, test='F', vcov=sandwich::vcovHC)
  # pbj methods
  statmap <- lmPBJ(pain$data$images, form = ~ ns(x, df = 4),
                   formred = ~ x, mask = mask,
                   template=pain$template, data = pain$data,
                   Winv = pain$data$Winv, zeros=TRUE)
  expect_equal(statmap$coef[,1], coefficients(model)[-c(1,2)], tolerance=tol )
  expect_equal(qf(pchisq(statmap$stat, df=statmap$df), df1 = statmap$df, df2=statmap$rdf)[1], waldtestres$F[2], tolerance=tol)
})

# check for errors
statmap <- lmPBJ(pain$data$images, form = ~ ns(x, df = 4),
                 formred = ~ x, mask = pain$mask,
                 template=pain$template, data = pain$data,
                 Winv = pain$data$Winv, zeros=TRUE)
pbjtest = pbjSEI(statmap, nboot = 5, cfts.s = c(0.1, 0.25))
pbjtest = pbjSEI(statmap, nboot = 5, cfts.p = c(0.01, 0.05))
