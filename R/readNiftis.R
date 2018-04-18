#' A Function to read NIfTI files in parallel
#'
#' This function is a wrapper for readNifti that opens a vector of nii or
#' nii.gz files in parallel.
#' @param files A vector of image names.
#' @param mc.cores number of cores to use when loading images
#' @keywords cats
#' @export
#' @examples
readNiftis = function(files=NULL, mc.cores=getOption("mc.cores", 2L)){
  if(length(files)==1){
    imgs = RNifti::readNifti(files)
  } else {
    imgs = parallel::mclapply(files, RNifti::readNifti, mc.cores=mc.cores)
    imgs = do.call(abind::abind, list(imgs, along=4))
  }
}
