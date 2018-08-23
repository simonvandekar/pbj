
sizing <- function(fn) paste0(round(file.info(fn)$size/1024^2, 1), "M")
statFile  <- function(label, fn) paste0(label, fn, "' ", sizing(fn), "\n")
statNifti <- function(label, im)
{
  paste0(
    label, "nifti[",
    paste0(dim(im), collapse=" x "),
    "] Pixel dimensions: ",
    paste0(attr(im, "pixdim"), c(attr(im, "pixunits")[1], attr(im, "pixunits")), collapse=" x "),
    '\n'
  )
}
statMatrix <- function(label, mm)
{
  paste0(
    label, "matrix[",
    paste0(dim(mm), collapse=" x "),
    "]\n"
  )
}

statInner <- function(label, obj)
{
  if(is.null(obj))    return("")
  if(all(is.na(obj))) return(paste0(label, "NA"))
  
  if(class(obj)[1] == "character")  return(statFile(label, obj))
  if(class(obj)[1] == "niftiImage") return(statNifti(label, obj))
  if(class(obj)[1] == "matrix")     return(statMatrix(label, obj))
  
  paste0(label, "Unhandled Class(",class(obj)[1],")")
}

#' @export
summary.statMap <- function(object, ...)
{
  paste0(
    "\nFormula: ", paste0(as.character(object$formulas[[2]]), collapse=''), paste0(as.character(object$formulas[[1]]), collapse=''), "\n",
    "\nFiles:\n",
    statInner("  Stat:       '", object$stat),
    statInner("  Sqrt Sigma: '", object$sqrtSigma),
    statInner("  Mask:       '", object$mask),
    statInner("  Template:   '", object$template),
    "  Robust:     ", object$robust
  )
}

#' @export
print.statMap <- function(x, ...)
{
  cat(summary(x, ...))
}

plot.statMap <- function(x, slice=1, ...)
{
  stop("FILL IN PLOTTING FUNCTION HERE. SLICE IS an example of an additional parameter")
}