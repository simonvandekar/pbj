# loading example data
library(pain21)
pain = pain21()
pdata = pain$data

# fitting intercept only model, weights proportional to study sample size
pbjModel1 = lmPBJ(images=pdata$images, form=~1, formred=~0, W = pdata$n,
                  mask=pain$mask, data=pdata, template = pain$template)
pbjModel1


# image(pbjModel1, index=5:15, nrow=2)
# image(pbjModel1)
# par(mfrow=c(1,2))
# image(pbjModel1, title='Pain Analysis', oma=c(0,0,3,0), crop=TRUE, index=2, lo=FALSE)
# image(pbjModel1, oma=c(0,0,3,0), crop=FALSE, index=2, lo=FALSE)
# image(stat.statMap(pbjModel1), title='Pain Analysis', oma=c(0,0,3,0), crop=FALSE, index=45)

# fitting regression of images onto study sample size, weights proportional to study sample size
pbjModel2 = lmPBJ(images=pdata$images, form=~n, formred=~1, W = pdata$n, mask=pain$mask, data=pdata)
pbjModel2
