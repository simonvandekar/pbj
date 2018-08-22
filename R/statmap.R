sizing <- function(fn) paste0(round(file.info(fn)$size/1024^2, 1), "M")

#' @export
summary.statMap <- function(object, ...)
{
  paste0(
    "\nFormula: ", paste0(as.character(object$formulas[[2]]), collapse=''), paste0(as.character(object$formulas[[1]]), collapse=''), "\n",
    "\nFiles:\n",
    "  Stat:       '", object$stat,      "' ", sizing(object$stat), "\n",
    "  Sqrt Sigma: '", object$sqrtSigma, "' ", sizing(object$sqrtSigma), "\n",
    "  Mask:       '", object$mask,      ", ", sizing(object$mask),"\n",
    if(is.null(object$template)) "" else paste0("  Template:  '", object$template,      "', ", sizing(object$template),"\n")
  )
}

#' @export
print.statMap <- function(x, ...)
{
  cat(summary(x, ...))
}