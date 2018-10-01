#' Builds (semi)Parametric Bootstrap Joint ((s)PBJ) Coverage Probability Excursion (CoPE) Sets
#'
#' @param statMap statMap object as obtained from computeStats.
#' @param rsq minimum effect size of interest. Currently, only supported for length(rsq)=1.
#' @param nboot Number of bootstraps to perform.
#' @param boundary Sample from the boundary as in Sommerfeld et al. (2018).
#'
#' @return Returns a list of length length(rsq)+2. The first two elements contain
#' statMap$stat and statMap$template. The remaining elements are lists containing the following:
#' \item{Aminus}{Aminus<=1-alpha contains E{statMap$stat}>chsq with probability 1-alpha.}
#' \item{Aplus}{Aplus>1-alpha is contained in E{statMap$stat}>chsq with probability 1-alpha.}
#' \item{rsq}{A minimum threshold to define interesting locations. It is the proportion of signal in the chi-square statistic.}
#' \item{chsq}{rsq converted to the chi-squared scale using chsq = ( (n-1) \* rsq + 1) /(1-rsq) \* df.}
#' \item{CDFs}{The bootstrap CDFs used to compute Aminus and Aplus.}
#' @export
#' @importFrom stats ecdf qchisq rnorm
#' @importFrom RNifti writeNifti updateNifti
pbjExSet = function(statMap, rsq=0.05, nboot=5000, boundary=FALSE, eps=0.01){
  if(class(statMap)[1] != 'statMap')
    warning('Class of first argument is not \'statMap\'.')
  
  mask = if(is.character(statMap$mask)) readNifti(statMap$mask) else statMap$mask
  stat = if(is.character(statMap$stat)) readNifti(statMap$stat) else statMap$stat
  template = statMap$template
  if(is.null(template)) template = mask
  df = statMap$df
  rdf = statMap$rdf
  
  if(df==0){
    sgnstat = sign(stat)
    stat = stat^2
    df=1; zerodf=TRUE
  } else {
    zerodf=FALSE
  }
  
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
  chsq = ( (n-1) * rsq + 1) /(1-rsq) * df
  # first column is min in stat>chsq
  # second column is max in stat<=chisq
  # boundary only uses the boundary voxels as in Sommerfeld et al. 2018
  if(boundary & df==1){
    bmask = which(stat[ mask!=0] <= chsq+eps & stat[ mask!=0]>= chsq-eps )
    stat = sqrt(stat) * sgnstat
    Aminus = Aplus = stat
    # only need the second column here
    # In this case set chisq=0, so that we are taking max over all voxels in the boundary
    Fs = pbjESboundary(sqrtSigma[bmask,], nboot)
    Fs = ecdf(Fs)
    Aminus[mask!=0] = 1-Fs( abs(Aminus[mask!=0]) - sqrt(chsq) )
    Aplus[mask!=0] = Fs(abs(Aplus[mask!=0]) - sqrt(chsq))
  } else {
    Aminus = Aplus = stat
    Fs = pbjES(stat[mask==1], sqrtSigma, chsq, df, nboot)
    Fs = apply(Fs, 2, ecdf)
    Aminus[mask!=0] = 1-Fs[[1]](Aminus[mask!=0])
    Aplus[mask!=0] = Fs[[2]](Aplus[mask!=0])
    if(zerodf){
      stat = sqrt(stat) * sgnstat
    }
  }
  
  out = list(stat=stat, template=template, rsq=list(Aminus=Aminus, Aplus=Aplus, rsq=rsq, chsq=chsq, CDFs=Fs) )
  names(out)[3] =  paste0('rsq', rsq)
  class(out) = c('CoPE', 'list')
  return(out)
}
