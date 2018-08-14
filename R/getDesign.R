#' Gets a Design Matrix for a GAM or Linear Model
#'
#' This function gets a design matrix for a generalized additive model (GAM)
#'  from the mgcv package or a linear model.
#' @param form A formula or string beginning with ~ that can be interpreted as an mgcv or lm
#'  formula without the outcome variable.
#' @param data A data frame containing the variables specified in form.
#' @keywords design matrix
#' @return Returns a design matrix constructed from form and data.
#' @importFrom mgcv gam
#' @importFrom stats as.formula model.matrix update.formula
#' @export
getDesign = function(form, data){
  # if it looks like a gam formula
  if(any(grepl("s\\(", form))){
    # 1:n because you need to give gam an outcome to easily get the design matrix
    data$x = 1:nrow(data)
    lmfull = if(class(form)=='formula') update.formula(form, x ~ .) else paste('x', form)
    lmfull = model.matrix(mgcv::gam(as.formula(lmfull), data=data) )

  # else it's a linear model
  } else {
    lmfull = model.matrix(as.formula(form), data=data)
  }
  return(lmfull)
}
