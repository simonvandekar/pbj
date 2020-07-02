#' Gets a Full and Reduced Design Matrices for a Linear Model
#'
#' This function gets full and reduced design matrix for a linear model as well as the degrees of freedom for the comparison.
#' @param form A formula or string beginning with ~ that can be interpreted as an lm
#'  formula without the outcome variable. Alternatively, a design matrix.
#' @param formred Similar to above. except this should be a simplified model relative to form.
#'  Put another way, the design matrix for this model should span a subset of the column space of the other model
#' @param data A data frame containing the variables specified in form.
#' @param tol Tolerance for determining the number of linearly independent columns to determine the df of the test.
#' @keywords design matrix
#' @return Returns design matrices for the full and reduced model and the df for the comparison between the two.
#' @details If robust=TRUE, then the number of parameters being tested has to
#' be equal to the df of the test, if not, then the covariance matrix of the parameters will be noninvertible.
#' If robust=TRUE, then the design is rotated to a df lower dimensional space. This problem happens when
#' testing splines, where the linear component of the full model is parameterized differently than the reduced model.
#' @importFrom stats as.formula model.matrix update.formula get_all_vars
#' @export
getDesign = function(form, formred, data, tol=1e-7){
  data=na.omit(get_all_vars(form, data=data))
  if(!is.matrix(form) & !is.matrix(formred)){
    X = model.matrix(as.formula(form), data)
    Xred = if(!is.null(formred)) model.matrix(as.formula(formred), data) else NULL
  } else {
    X = form
    Xred = formred
    form <- formred <- NULL
  }

  X.svd = svd(qr.resid(qr(Xred), X), nv=0 )
  # using svd to get full model df. Accounts for the possibility that some columns of Xred are linearly dependent on X.
  # For example, with ns using spline basis functions.
  df = sum(X.svd$d/sum(X.svd$d)>tol)
  cols = sum(!colnames(X) %in% colnames(Xred))
  if(df< cols ){
    message('df=',df, ' is less than additional number of columns in full model (', cols,
            '). \nCoefficients will likely be uninterpretable.' )
      message('Creating new lower dimensional basis with df=', df, '.')
      X1 = X.svd$u[,1:df]
      colnames(X1) = paste0('u', 1:df)
      X = cbind(Xred, X1)
  }
  return(list(X=X, Xred=Xred, df=df))
}
