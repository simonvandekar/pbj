# gets pvalues for null and alternative clusters
getPvalues = function(pmap=NULL, betamap=NULL, kernel='box', pfunc=function(x) 10^(-x) ){
	if(is.null(pmap))
		stop("pmap is required.")
	if(is.character(pmap))
		pmap = readNifti(pmap)
	if(is.character(betamap))
		betamap = readNifti(betamap)

	# get betamap components
	betacomp = mmand::components(sig, mmand::shapeKernel(3, 3, type='box'))
	clustsize = sort(table(c(betacomp)))
	names(clustsize) = 1:length(clustsize)

	# renumber components smallest to largest
	for(i in 1:length(clustsize)){
		betacomp[ which(betacomp==as.numeric(names(clustsize))[i]) ] = i + length(clustsize)
	}
	betacomp = betacomp - length(clustsize)

	comps = mmand::components(pmap, mmand::shapeKernel(3, 3, type=kernel))
	if(!is.null(betacomp)){
		allaltinds = list()
		altp = list()
		# compute p-values for each real cluster defined by betacomp
		for(clind in 1:max(betacomp, na.rm=TRUE)){
			altinds = na.omit(unique(comps[ betacomp==clind ]))
			altp[[clind]] = if(length(altinds)==0) 1 else sort(pfunc(sapply(altinds, function(x) pmap[ which(comps==x)[1] ] )) )
			allaltinds[[clind]] = altinds
		}
	} else {
		altp <- altind <- NULL
	}
	nullinds = na.omit(unique(c(comps)))
	nullinds = nullinds[ !nullinds %in% unlist(allaltinds) ]
	nullp = if(length(nullinds)==0) 1 else sort(pfunc(sapply(nullinds, function(x) pmap[ which(comps==x)[1] ] )) )
	list(altp = altp, clustsize = clustsize, nullp = nullp, altinds = allaltinds, componentmap=comps)
}
