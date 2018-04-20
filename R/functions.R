### RANDOMISE WRAPPER ###
randomise = function(filelist=NULL, dat=NULL, maskfile=NULL, X=NULL, Xred=NULL, outdir=NULL, groupvar=NULL, nperm=500, thresh=0.01, run=TRUE){
	outdirdir = dirname(outdir)
	mergednifti=file.path(outdirdir, 'merged.nii.gz')
	if(length(filelist)==1){
		mergednifti=filelist
	}
	# merge images
	if(!file.exists(mergednifti)){
		mergeimages(filelist, mergednifti)
	}

	## DESIGN AND CONTRASTS ##
	# design file
	n = nrow(X)
	p = ncol(X)
	# If there is only 1 column then do a one-sample t-test
	if(p>1){
		p2 = p - ncol(Xred)
		matfile = paste(outdir, 'design.mat', sep='_')
		cat('/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(X), '\n/PPheights\t', paste(apply(X, 2, function(x) abs(diff(range(x))) ), collapse='\t'), '\n\n/Matrix\n', sep='', file=matfile)
		write.table(X, append=TRUE, file=matfile, row.names=FALSE, col.names=FALSE)

		# contrast file
		confile1 = paste(outdir, 'design.con', sep='_') # for f-test
		cons = matrix(0, nrow=p2, ncol=ncol(X))
		cons[ cbind(1:(p2), which(! colnames(X) %in% colnames(Xred) ) ) ] = 1
		cat('/ContrastName1\t temp\n/ContrastName2\t\n/NumWaves\t', ncol(X), '\n/NumPoints\t', nrow(cons), '\n/PPheights\t', paste(rep(1,ncol(cons)), collapse='\t'), '\n/RequiredEffect\t1\t1\n\n/Matrix\n', sep='', file=confile1)
		write.table(cons, append=TRUE, file=confile1, row.names=FALSE, col.names=FALSE)

		# fts file
		ftsfile = paste(outdir, 'design.fts', sep='_')
		fts = matrix(1, nrow=1, ncol=nrow(cons)) # ftest of all contrasts
		cat('/NumWaves\t', nrow(cons), '\n/NumContrasts\t', 1, '\n\n/Matrix\n', sep='', file=ftsfile)
		write.table(fts, append=TRUE, file=ftsfile, row.names=FALSE, col.names=FALSE)

		fcmd = paste('randomise -i', mergednifti, '-o', outdir, '-d', matfile, '-t', confile1, '-f', ftsfile, '--fonly -F', qf( (1-thresh),df1=p2, df2=(n-p) ), '-x -N -T -n', nperm, '--uncorrp' )

		# grp file UNTESTED
		if(!is.null(groupvar)){
			grpfile = paste(outdir, "design.grp", sep='_')
			cat('/NumWaves\t1\n/NumPoints\t', nrow(dat), '\n/Matrix\n', file=grpfile)
			write.table(groupvar, append=TRUE, file=grpfile, row.names=FALSE, col.names=FALSE)
			fcmd = paste(fcmd, '-e', grpfile)
		}

		if(!is.null(maskfile)){
			fcmd = paste(fcmd, '-m', maskfile)
		}

	} else {
		fcmd = paste('randomise -i', mergednifti, '-o', outdir, '-1 -c', qt( (1-thresh),df=(n-1) ), '-x -N -T -n', nperm, '--uncorrp' )
		if(!is.null(maskfile)){
			fcmd = paste(fcmd, '-m', maskfile)
		}
	}
		# randomise saves everything-- just return the commands
		if(run){
			system(fcmd)
		}
		return(fcmd)
}

# merges a list of images using fslmerge
mergeimages = function(filelist=NULL, mergednifti=NULL){
		txt = tempfile( tmpdir="./")
		cat('merging images\n')
		cmd = paste('fslmerge -t', mergednifti, paste(filelist, collapse=' '))
		cat(cmd, file=txt)
		# for some reason I can't apply system to this command directly
		# I think due to limitations on the string length system can take
		# Instead, I cat to a txt file, and then use system to evaluate the
		# text file
		mergt = system.time(trash <- system(paste("$(cat", txt, ")") ) )[3]
		unlink(txt)
		cat('merge time:', (mergt)/60, 'minutes\n')
}


# specify outroot to keep output
easythresh = function(pstat=NULL, mask=NULL, cthresh=0.01, pthresh=0.15, bgimg=mask, run=TRUE, outroot=NULL){
	if(any(is.null(c(pstat, mask)))){
		stop('Must specify pstat and mask')
	}
	if(is.null(outroot)){
		outroot0 = file.path(tempfile(), 'et')
	} else {outroot0=outroot}
	dir.create(dirname(outroot0), showWarnings=FALSE)
	zstat = file.path(dirname(outroot0), 'zstat.nii.gz')
	system(paste('fslmaths', pstat, '-sub 1 -mul -1 -ptoz', zstat))
	here = getwd()
	etcmd = paste('easythresh', zstat, mask, qnorm(cthresh, lower.tail=FALSE), pthresh, bgimg, basename(outroot0))
	#
	if(run){
	setwd(dirname(outroot0))
		system(etcmd)
		ptab = read.table(paste('cluster_',basename(outroot0), '.txt', sep=''), header=TRUE, sep='\t')
		if(nrow(ptab)==0){
			pvalues=NA
			clustermask=NA
		} else {
			pvalues = ptab[,'P']
			clustermask = readNifti(paste('cluster_mask_', basename(outroot0), '.nii.gz', sep='') )
			pmap = clustermask
			Max = max(clustermask)
			for(clind in 1:Max){
				pmap[ pmap==(Max - clind+1) ] = pvalues[ clind]
			}
		}
	setwd(here)
	}
	if(is.null(outroot)){
		unlink(dirname(outroot0), recursive=TRUE)
	}
	out = list(pvalues=pvalues, cmd = etcmd, clustermask=clustermask, pmap=pmap)
}

regionmean = function(input=NULL, label=NULL){
	vec = system(paste('c3d', input, label, '-lstat'), intern=TRUE)
	# remove zero label get first two columns
	vec = read.table(text=vec, strip.white=TRUE, skip=1)[-1,1:2]
}

