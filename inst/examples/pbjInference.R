# loading example data
library(pain21)
pain = pain21()
pdata = pain$data

# fitting regression of images onto study sample size, weights proportional to study sample size
pbjModel2 = lmPBJ(images=pdata$images, form=~n, formred=~1, W = pdata$n, mask=pain$mask, data=pdata)
pbjModel2

cft = qchisq(0.01, df=pbjModel2$sqrtSigma$df, lower.tail=FALSE)
pbjModel2 = pbjInference(pbjModel2, nboot=10, cft=cft)
pbjModel2

