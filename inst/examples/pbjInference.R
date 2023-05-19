# loading example data
library(pain21)
pain = pain21()
pdata = pain$data

# fitting regression of images onto study sample size, weights proportional to study sample size
pbjModel2 = lmPBJ(images=pdata$images, form=~n, formred=~1, W = pdata$n, mask=pain$mask, data=pdata)
table.statMap(pbjModel2, cft_p=0.01)


# p-value thresholding for cluster extent/mass
# run in parallel
pbjModelAll <- pbjInference(pbjModel2, nboot=50, cft_p=0.01, CEI=TRUE, CMI=TRUE, max=TRUE, mc.cores=ceiling(parallel::detectCores()/2))
head(table.statMap(pbjModelAll, method = 'maxima'))
head(table.statMap(pbjModelAll, method = 'CEI', cft_p=0.01))
# returns the first threshold used
head(table.statMap(pbjModelAll, method = 'CMI'))
# RESI effect size thresholding for cluster extent/mass
pbjModel2 = pbjInference(pbjModel2, nboot=2, cft_s=c(0.1, 0.25), CMI=TRUE)
head(table.statMap(pbjModel2, method = 'CEI', cft_s=0.1))
head(table.statMap(pbjModel2, method = 'CEI', cft_s=0.25))
# Inference on local maxima
pbjModel2 = pbjInference(pbjModel2, nboot=2, statistic = maxima)
# Cluster extent inference
pbjModelCEI = pbjInference(pbjModel2, nboot=2, statistic = cluster, cft_s=0.1)

# not run
outfile = paste0(tempfile(), '.rds')
bgInfo <- pbjInference(pbjModel2, rdata_rds = outfile, nboot=50, cft_p=0.01, CEI=TRUE, CMI=TRUE, max=FALSE, mc.cores=ceiling(parallel::detectCores()/2))
