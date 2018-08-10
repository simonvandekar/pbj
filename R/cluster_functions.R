
# cluster correction
# c.values is a vector of cluster extents
# residuals are nXV matrix where V is number of voxels in the image
# thr is a vector of cluster forming thresholds as a p-value
# kernel is for connected components definition. "box" or "diamond" are both reasonable

#' @importFrom stats quantile rnorm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom RNifti readNifti
#' @importFrom mmand shapeKernel
pbjES = function(obs, residuals=NULL, mask, df=1, rdf, alpha=0.05, thr=df*rdf/(rdf-2), nsim=5000){

	if(is.null(thr))	     stop('Must specify cluster forming threshold.')
	if(is.character(mask)) mask = readNifti(mask)
	if(is.character(obs))  obs  = readNifti(obs)

	k = shapeKernel(3, 3, type='diamond')
	thrmask = mask
	thrmask = ( obs > thr)
	# gets boundary voxels
	thrmask = dilate(thrmask, k) - thrmask
	# get index vector for residuals
	thrmask = thrmask[ mask==1]

	if(is.character(residuals)){
		residuals = readNifti(residuals)
	}
	residuals = apply(residuals, 4, function(x) x[mask==1])

	ssqs = sqrt(rowSums(residuals^2))
	if( any(ssqs != 1) ){
		residuals = sweep(residuals, 1, ssqs, '/')
	}
	p=nrow(residuals)
	r = min(rdf, p)
	residuals = svd(residuals, nu=r, nv=0)
	residuals = sweep(residuals$u, 2, residuals$d[1:r], "*")
	# get residuals only on the boundary
	residuals = residuals[ thrmask,]

	Fs = matrix(NA, nsim, 2)
	pb = txtProgressBar(style=3, title='Generating null distribution')
	for(i in 1:nsim){
		S = matrix(rnorm(r*df), r, df)
		statimg = rowSums((residuals %*% S)^2)
		Fs[i,] = c(min(statimg), max(statimg))
		setTxtProgressBar(pb, round(i/nsim,2))
	}
	close(pb)
	# this set is contained in the actual set with probability 1-\alpha/2
	ESsub = obs * ((obs - thr) > (quantile(Fs[,2], 1-alpha/2) - df) )
	# this set contains the actual set with probability 1-alpha/2
	ESsup = obs * ( (obs - thr) > (quantile(Fs[,1], alpha/2) - df) )
	return(list(ESsub, ESsup))
}
