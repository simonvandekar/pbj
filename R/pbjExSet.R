#' Constructs (semi)Parametric Bootstrap Joint ((s)PBJ) Coverage Probability Excursion (CoPE) Sets
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param ses minimum effect size of interest for.
#' @param nboot Number of bootstraps to perform.
#' @param boundary Sample from the boundary as in Sommerfeld et al. (2018).
#' @param eps tolerance for the chi-square statistic.  Only used if
#' \code{boundary} is \code{TRUE} and the degrees of freedom is 1.
#'
#' @return Returns a list of length length(ses)+2. The first two elements contain
#' statMap$stat and statMap$template. The remaining elements are lists containing the following:
#' \item{Aminus}{Aminus<=1-alpha contains E{statMap$stat}>chsq with probability 1-alpha.}
#' \item{Aplus}{Aplus>1-alpha is contained in E{statMap$stat}>chsq with probability 1-alpha.}
#' \item{ses}{A minimum threshold to define interesting locations. It is the proportion of signal in the chi-square statistic.}
#' \item{chsq}{ses converted to the chi-squared scale using chsq = ( (n-1) \* ses + 1) /(1-ses) \* df.}
#' \item{CDFs}{The bootstrap CDFs used to compute Aminus and Aplus.}
#' @export
#' @importFrom stats ecdf qchisq rnorm
#' @importFrom RNifti writeNifti updateNifti
pbjExSet = function(statMap, ses=0.2, nboot=5000, boundary=FALSE, eps=0.01){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')

  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  stat = if(is.character(statMap$stat)) readNifti(statMap$stat)[ mask!=0] else statMap$stat
  template = statMap$template
  if(is.null(template)) template = mask
  df = statMap$df
  rdf = statMap$rdf

  # load sqrtSigma if nifti image
  if(is.character(statMap$sqrtSigma)){
    sqrtSigma = readNifti(statMap$sqrtSigma)
    sqrtSigma = apply(sqrtSigma, 4, function(x) x[mask==1])
  } else {
    sqrtSigma = statMap$sqrtSigma
    rm(statMap)
  }

  # normalize sqrtSigma
  ssqs = sqrt(rowSums(sqrtSigma^2))
  if( any(ssqs != 1) ){
    sqrtSigma = sweep(sqrtSigma, 1, ssqs, '/')
  }
  n = ncol(sqrtSigma)
  chsq_threshold = n * ses + df
  Aminus = Aplus = mask
  # first column is min in stat>chsq
  # second column is max in stat<=chisq
  # boundary only uses the boundary voxels as in Sommerfeld et al. 2018
  if(boundary & df==0){
    bmask = which(stat <= sqrt(chsq_threshold)+eps & stat >= sqrt(chsq_threshold)-eps )
    # only need the second column here
    # In this case set chisq=0, so that we are taking max over all voxels in the boundary
    Fs = pbjESboundary(sqrtSigma[bmask,], nboot)
    Fs = ecdf(Fs)
    Aminus[mask!=0] = 1-Fs( Aminus[mask!=0 ] + sqrt(chsq_threshold-df) )
    Aplus[mask!=0] = 1-Fs( Aplus[mask!=0] - sqrt(chsq_threshold-df))
  } else if(!boundary & df==0){
    Fs = pbjESzerodf(stat[mask!=0], sqrtSigma, sqrt(chsq_threshold - df), nboot)
    Fs = apply(Fs, 2, ecdf)
    Aminus[mask!=0] = 1-Fs[[1]](Aminus[mask!=0])
    Aplus[mask!=0] = Fs[[2]](Aplus[mask!=0])
  } else {
    Aminus = Aplus = mask
    Fs = pbjES(stat, sqrtSigma, chsq_threshold, df, nboot)
    Fs = apply(Fs, 2, ecdf)
    Aminus[mask!=0] = 1-Fs[[1]](stat)
    Aplus[mask!=0] = Fs[[2]](stat)
    if(zerodf){
      stat = sqrt(stat) * sgnstat
    }
  }

  out = list(stat=stat, template=template, ses=list(Aminus=Aminus, Aplus=Aplus, ses=ses, CDFs=Fs) )
  names(out)[3] =  paste0('ses', ses)
  class(out) = c('CoPE', 'list')
  return(out)
}
