#' Gets p-values for Simulation Study
#'
#' @param pmap A cluster-extent p-value map (nii, nii.gz, or Nifti object),
#'  where each cluster is labeled with its FWE adjusted cluster-extent p-value.
#'  (Required argument)
#' @param betamap A known betamap (nii, nii.gz, or Nifti object) used in simulations to generate artificial
#'  signal in the imaging outcome. If NULL then power statistics are not
#'  reported.
#' @param kernel What type of kernel to use for connected components. Defaults
#'  to "box".
#' @param pfunc If pmap is transformed p-values then what function should be
#'  used to compute untransformed p-values. Uses pbj default, which stores p-
#'  values as -log10(p).
#'
#' @return Returns a list with the following values:
#' \item{altp}{If betamap is non-NULL then this returns a list the length of
#'  the number of connected components in betamap. Each list element gives
#'  the ordered p-values of clusters in pmap that have nonempty overlap with
#'  that cluster in betamap. The output is ordered by increasing size of
#'  connected components in betamap. NULL if betamap is NULL.}
#' \item{nullp}{If betamap is non-NULL then these are the p-values for clusters
#'  in pmap that have empty intersection with betamap. Otherwise, this is a
#'  vector of p-values for all clusters in pmap.}
#' \item{clustsize}{If betamap us non-NULL then this is the sizes of clusters in betamap in increasing order. NULL otherwise.}
#' \item{altinds}{If betamap is non-NULL then this returns a list the length of
#'  the number of connected components in betamap. Each list element gives the
#'  indices of clusters indices in componentmap that overlap with the
#'  corresponding cluster in betamap. NULL otherwise.}
#' \item{nullinds}{Vector of indices of components in component map that do not
#'  overlap with betamap.}
#' \item{componentmap}{Nifti object giving connected components of pmap.}
#' @export
#' @importFrom stats na.omit
#' @importFrom RNifti readNifti
#' @importFrom mmand shapeKernel
getPvalues = function(pmap=NULL, betamap=NULL, kernel='box', pfunc=function(x) 10^(-x) ){
	if(is.null(pmap))
		stop("pmap is required.")
	if(is.character(pmap))
		pmap = readNifti(pmap)
	if(is.character(betamap))
		betamap = readNifti(betamap)

	# get pmap components
	comps = mmand::components(pmap, mmand::shapeKernel(3, 3, type=kernel))

	if(!is.null(betamap)){
		# get betamap components
		betacomp = mmand::components(betamap, mmand::shapeKernel(3, 3, type=kernel))
		clustsize = sort(table(c(betacomp)))

		# renumber components smallest to largest
		for(i in 1:length(clustsize)){
			betacomp[ which(betacomp==as.numeric(names(clustsize))[i]) ] = i + length(clustsize)
		}
		names(clustsize) = 1:length(clustsize)
		betacomp = betacomp - length(clustsize)

		allaltinds = list()
		altp = list()
		# compute p-values for each real cluster defined by betacomp
		for(clind in 1:max(betacomp, na.rm=TRUE)){
			altinds = na.omit(unique(comps[ betacomp==clind ]))
			altp[[clind]] = if(length(altinds)==0) 1 else sort(pfunc(sapply(altinds, function(x) pmap[ which(comps==x)[1] ] )) )
			allaltinds[[clind]] = altinds
		}
	} else {
		clustsize <- altp <- allaltinds <- NULL
	}
	nullinds = na.omit(unique(c(comps)))
	nullinds = nullinds[ !nullinds %in% unlist(allaltinds) ]
	nullp = if(length(nullinds)==0) 1 else sort(pfunc(sapply(nullinds, function(x) pmap[ which(comps==x)[1] ] )) )
	list(altp = altp, clustsize = clustsize,  altinds = allaltinds, nullp = nullp, nullinds=nullinds, componentmap=updateNifti(comps, pmap))
}
