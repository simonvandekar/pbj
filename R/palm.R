#' this is an R wrapper function for PALM
#' it has limited functionality for the features that PALM offers
#'
#' @param i vector of strings pointing to input images or a single 4d nifti string.
#' @param m mask file.
#' @param s vector of strings pointing to input surfaces or a single 4d surface.
#' @param d design matrix as a data frame or csv file.
#' @param o output directory and prefix.
#' @param C cluster forming threshold as a Z-value
#' @param t T-contrast matrix that randomise takes as matrix, data frame, or csv.
#' @param f F-contrast matrix that randomise takes, as matrix data frame, or csv.
#' @param fonly only run the F-tests.
#' @param n number of permutations to run.
#' @param eb exchangeability blocks (not functional).
#' @param within not functional.
#' @param whole not functional.
#' @param ee exchangeable errors assumed.
#' @param T Run TFCE.
#'
#' @return Returns a string of the command that was evaluated.
#' @export
palm = function(i,m,s,d, o, C, t=NULL,f=NULL,fonly=TRUE, n=10000, eb=NULL, within=NULL, whole=NULL, ee=TRUE, T=FALSE){
  # if multiple arguments are passed assume it is a list of 3d niftis
  if(length(i)!=1){
     temp = paste0(tempfile(), '.nii.gz')
     mergeimages(i, temp)
     i = temp
  }
  # convert data frame to csv file
  if(is.data.frame(d)){
     temp = paste0(tempfile(), '.csv')
     write.table(d, temp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
     d = temp
  }
  X = read.csv(d,)

  # T-contrasts
  if(is.null(t)){
    t = diag(ncol(X))
  }
  if(is.matrix(t) | is.data.frame(t)){
    temp = paste0(tempfile(), '.csv')
    write.table(as.matrix(t), temp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    t = temp
  }

  # F-contrasts
  if(is.null(f)){
    f = matrix(c(0, rep(1,ncol(X)-1)), nrow=1)
  }
  if(is.matrix(f) | is.data.frame(f)){
    temp = paste0(tempfile(), '.csv')
    write.table(as.matrix(t), temp, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    f = temp
  }

  cmd = 'palm'
  argnames = names(match.call())[-1]
  args = lapply(argnames, function(x) eval(parse(text=x)))
  names(args) = argnames
  for( arg in names(args)){
    if(!is.null(args[[arg]])){
      if(is.logical(args[[arg]])){
        cmd = paste0(cmd, ' -', arg)
      } else {
        cmd = paste0(cmd, ' -', arg, ' ', args[[arg]])
      }
    }
  }
  system(cmd)
  cmd
}
