#' Gets a Design Matrix for a GAM or Linear Model
#'
#' This function gets a design matrix for a generalized additive model (GAM)
#'  from the mgcv package or a linear model.
#' @param eq A formula or string beginning with ~ that can be interpreted as an mgcv or lm
#'  formula without the outcome variable.
#' @param data A data frame containing the variables specified in eq.
#' @keywords design matrix
#' @return Returns a design matrix constructed from eq and data.
#' @export
getDesign = function(eq=NULL, data=NULL){
  # if it looks like a gam formula
  if(any(grepl("s\\(", eq))){
    # 1:n because you need to give gam an outcome to easily get the design matrix
    x = 1:nrow(data)
    lmfull = if(class(eq)=='formula') update.formula(eq, x ~ .) else paste('x', eq)
    lmfull = model.matrix(mgcv::gam(as.formula(lmfull), data=data) )

  # else it's a linear model
  } else {
    lmfull = model.matrix(as.formula(eq), data=data)
  }
  return(lmfull)
}
