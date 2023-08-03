#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>

static SEXP pbj_useLAPACK_Sym = NULL;

/* See https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Handling-lists */
static SEXP pbj_getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

SEXP pbj_pbjBootRobustX(SEXP qr, SEXP res, SEXP x1res, SEXP idmat, SEXP h, SEXP df) {
  /* Reproduce this:

    rsd <- qr.resid(sqrtSigma$QR, sqrtSigma$res)
    rsd <- sweep(rsd, 1, 1-h, '/')
    baz <- rep(list(rsd), df)
    qux <- simplify2array(baz)
    corge <- sweep(qux, c(1,3), sqrtSigma$X1res, '*')
    grault <- function(x){
      a <- qr.R(qr(x))
      b <- diag(ncol(x))
      backsolve(r=a, x=b)
    }
    garply <- apply(corge, 2, grault)
    BsqrtInv <- matrix(garply, nrow=df^2, ncol=V)

    foo <- function(ind) {
      bar <- crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind])
      baz <- matrix(BsqrtInv[,ind], nrow=df, ncol=df)
      crossprod(baz, bar)
    }
    qux <- lapply(1:V, foo)
    corge <- simplify2array(qux, higher=TRUE)
    statimg <- matrix(corge, nrow=df, ncol=V)
  */

  SEXP qr_qr, qr_qraux, elt, attr, dim, statimg;
  int *dim_ii, arr_idx_i, idval_i, k_i, n_i, ny_i, nrow_i, ncol_i, idncol_i,*idmat_ii, row_i, col_i, row2_i, col2_i,
      rsd_idx_i, idmat_idx_i, x1res_idx_i, corge_idx_i, df_i, layer_i, dim_prod_i, x_idx_i,
      ldx_i, p_i, *pivot_ii, bsqrtinv_idx_i, idx_i, df_sq_i, one_i;
  double *rsd_dd, *h_dd, *x1res_dd, *corge_dd, *idres_dd, *res_dd, *res2_dd, *x_dd, tol_d,
         *qraux_dd, *work_dd, *a_dd, *bsqrtinv_dd, one_d, *statimg_dd, zero_d;

  /* Type checking for qr */
  if (!isNewList(qr) || !inherits(qr, "qr")) {
    error("qr must be a QR decomposition");
  }

  if (isComplex(qr)) {
    error("not implemented for complex qr");
  }

  if (pbj_useLAPACK_Sym == NULL) {
    pbj_useLAPACK_Sym = install("useLAPACK");
  }

  attr = getAttrib(qr, pbj_useLAPACK_Sym);
  if (isLogical(attr) && length(attr) == 1 && LOGICAL(attr)[0] == TRUE) {
    error("not supported for LAPACK QR");
  }

  elt = pbj_getListElement(qr, "rank");
  if (isInteger(elt) && length(elt) == 1) {
    k_i = INTEGER(elt)[0];
  } else {
    error("qr$rank must be an integer vector of length 1");
  }

  qr_qr = pbj_getListElement(qr, "qr");
  if (isReal(qr_qr)) {
    dim = getAttrib(qr_qr, R_DimSymbol);
    if (isInteger(dim) && length(dim) == 2) {
      /* nrow(qr$qr) */
      n_i = INTEGER(dim)[0];
    } else {
      error("qr$qr must be a real 2d matrix");
    }
  } else {
    error("qr$qr must be a real 2d matrix");
  }

  qr_qraux = pbj_getListElement(qr, "qraux");
  if (!isReal(qr_qraux)) {
    error("qr$qraux must be a real vector");
  }

  /* Type checking for res */
  if (isReal(res)) {
    dim = getAttrib(res, R_DimSymbol);
    if (isInteger(dim) && length(dim) == 2) {
      dim_ii = INTEGER(dim);

      /* check nrow(res) == nrow(qr$qr) */
      nrow_i = dim_ii[0];
      if (nrow_i != n_i) {
        error("qr and res must have the same number of rows");
      }

      /* ncol(res) */
      ncol_i = ny_i = dim_ii[1];
    } else {
      error("res must be a real 2d matrix");
    }
  } else {
    error("res must be a real 2d matrix");
  }

  /* Type checking for idmat */
  /* Type checking for h */
  if(idmat != R_NilValue){
    if (!isInteger(idmat) || length(idmat) != n_i) {
      error("id must be NULL or an integer vector with the same length as nrow(res)");
    }
  }

  /* Type checking for df */
  if (!isInteger(df) || length(df) != 1) {
    error("df must be an integer vector of length 1");
  }
  df_i = INTEGER(df)[0];
  df_sq_i = df_i * df_i;

  /* Type checking for x1res */
  if (isReal(x1res)) {
    dim = getAttrib(x1res, R_DimSymbol);
    if (isInteger(dim) && length(dim) == 2) {
      dim_ii = INTEGER(dim);

      /* check nrow(x1res) == nrow(res) */
      if (dim_ii[0] != nrow_i) {
        error("res and x1res must have the same number of rows");
      }

      /* check ncol(x1res) == df */
      if (dim_ii[1] != df_i) {
        error("nrow(x1res) must equal df");
      }
    } else {
      error("x1res must be a real 2d matrix");
    }
  } else {
    error("x1res must be a real 2d matrix");
  }

  /* Type checking for h */
  if (!isReal(h) || length(h) != n_i) {
    error("h must be a real vector with the same length as nrow(res)");
  }

  /* Sanity check for qr.default/dqrdc2 */
  if (nrow_i * df_i > 2147483647) {
    error("too large a matrix for LINPACK");
  }

  /* Sanity check for convenience, not technical reasons */
  if (df_i > nrow_i) {
    error("df cannot be larger than nrow(res)");
  }

  /* Allocate space for rsd vector */
  rsd_dd = Calloc(length(res), double);

  /* Make a copy of res, since dqrrsd changes some values in place */
  res_dd = REAL(res);
  res2_dd = Calloc(length(res), double);
  memcpy(res2_dd, res_dd, length(res) * sizeof(double));

  /* rsd <- qr.resid(sqrtSigma$QR, sqrtSigma$res) */
  F77_CALL(dqrrsd)(REAL(qr_qr), &n_i, &k_i, REAL(qr_qraux), res2_dd, &ny_i, rsd_dd);
  Free(res2_dd);

  /* rsd <- sweep(rsd, 1, 1-h, '/') */
  h_dd = REAL(h);
  rsd_idx_i = 0;
  for (col_i = 0; col_i < ncol_i; col_i++) {
    for (row_i = 0; row_i < nrow_i; row_i++) {
      rsd_dd[rsd_idx_i] = rsd_dd[rsd_idx_i] / (1 - h_dd[row_i]);
      rsd_idx_i++;
    }
  }

  /* Allocate space for corge vector */
  corge_dd = Calloc(length(res) * df_i, double);

  /* baz <- rep(list(rsd), df)
     qux <- simplify2array(baz)
     corge <- sweep(qux, c(1,3), sqrtSigma$X1res, '*') */
  x1res_dd = REAL(x1res);
  x1res_idx_i = 0;
  corge_idx_i = 0;
  for (layer_i = 0; layer_i < df_i; layer_i++) {
    rsd_idx_i = 0;

    for (col_i = 0; col_i < ncol_i; ) {
      for (row_i = 0; row_i < nrow_i; row_i++) {
        corge_dd[corge_idx_i] = rsd_dd[rsd_idx_i] * x1res_dd[x1res_idx_i];
        corge_idx_i++;
        rsd_idx_i++;
        x1res_idx_i++;
      }

      col_i++;

      /* Reset x1res index back to the first row in its current column, unless
       * we're about to break out of this inner loop */
      if (col_i < ncol_i) {
        x1res_idx_i -= nrow_i;
      }
    }
  }
  Free(rsd_dd);



  /* Check if idmat is null if not multiply with corge_dd */
  if(idmat != R_NilValue){

    /* This creates a copy of idmat */
    idmat_ii = INTEGER(idmat);
    /* This will be nrows of res */

    /* get the maximum value of idmat */
    idncol_i = 0;
    for(idmat_idx_i=0; idmat_idx_i < nrow_i; idmat_idx_i++){
      idval_i = idmat_ii[idmat_idx_i];
      if(idncol_i < idval_i){
        idncol_i = idval_i;
      }
    }

    /* Allocate idres_dd. idncol_i X V X df array*/
    idres_dd = Calloc(idncol_i * ncol_i * df_i, double);

    idmat_idx_i = 0;
    corge_idx_i = 0;
    arr_idx_i=0;
    for (layer_i = 0; layer_i < df_i; layer_i++) {
      for (col_i = 0; col_i < ncol_i; col_i++) {
        for (row_i = 0; row_i < nrow_i; row_i++) {
          idmat_idx_i = (idmat_ii[row_i]-1) + arr_idx_i;
          idres_dd[idmat_idx_i] = idres_dd[idmat_idx_i] + corge_dd[corge_idx_i];
          corge_idx_i++;
        }
        /* idncol_i is number of rows of result of idres_dd */
        /* jump to next column*/
        arr_idx_i = arr_idx_i + idncol_i;
      }
    }



    /* replace corge_dd with idres_dd */
    nrow_i = idncol_i;
    Free(corge_dd);
    corge_dd = idres_dd;
    /*
    printf("%u\n", nrow_i);
    for(arr_idx_i=0; arr_idx_i< df_i*idncol_i*ncol_i; arr_idx_i++){
      printf("%f", corge_dd[arr_idx_i]);
      printf(" ");
    }
     */
  }

  /*
    grault <- function(x){
      a <- qr.R(qr(x))
      b <- diag(ncol(x))
      backsolve(r=a, x=b)
    }
    garply <- apply(corge, 2, grault)
  */
  x_dd = Calloc(nrow_i * df_i, double);
  qraux_dd = Calloc(df_i, double);
  pivot_ii = Calloc(df_i, int);
  work_dd = Calloc(df_i * 2, double);
  a_dd = Calloc(df_sq_i, double);
  bsqrtinv_dd = Calloc(df_sq_i * ncol_i, double);

  dim_prod_i = ncol_i * nrow_i;
  bsqrtinv_idx_i = 0;
  for (col_i = 0; col_i < ncol_i; col_i++) {
    /* create subset 2d matrix of corge */
    x_idx_i = 0;
    for (layer_i = 0; layer_i < df_i; layer_i++) {
      memcpy(x_dd + x_idx_i, corge_dd + (layer_i * dim_prod_i) + (col_i * nrow_i), nrow_i * sizeof(double));
      x_idx_i += nrow_i;
    }

    /* x <- qr(x) */
     /* SNV edits: used to set n_i=nrow_i. nrow_i is now equal to idncol_i. n_i is nrow of res */
    ldx_i = nrow_i;
    /* n_i = nrow_i; */
    p_i = df_i;
    tol_d = 0.0000001;
    k_i = 0;
    for (idx_i = 0; idx_i < df_i; idx_i++) qraux_dd[idx_i] = 0;
    for (idx_i = 0; idx_i < df_i; idx_i++) pivot_ii[idx_i] = idx_i + 1;
    for (idx_i = 0; idx_i < df_i * 2; idx_i++) work_dd[idx_i] = 0;
    F77_CALL(dqrdc2)(x_dd, &ldx_i, &nrow_i, &p_i, &tol_d, &k_i, qraux_dd, pivot_ii, work_dd);

    /* a <- qr.R(x) */
    idx_i = 0;
    for (col2_i = 0; col2_i < df_i; col2_i++) {
      x_idx_i = col2_i * nrow_i;
      for (row2_i = 0; row2_i < df_i; row2_i++) {
        if (row2_i > col2_i) {
          a_dd[idx_i] = 0;
        } else {
          a_dd[idx_i] = x_dd[x_idx_i];
        }
        idx_i++;
        x_idx_i++;
      }
    }

    /* Put an identity matrix into bsqrtinv_dd at the proper place */
    idx_i = bsqrtinv_idx_i;
    for (col2_i = 0; col2_i < df_i; col2_i++) {
      for (row2_i = 0; row2_i < df_i; row2_i++) {
        if (col2_i == row2_i) {
          bsqrtinv_dd[idx_i] = 1;
        } else {
          bsqrtinv_dd[idx_i] = 0;
        }
        idx_i++;
      }
    }

    /* bs <- backsolve(r=a, x=b) */
    one_d = 1.0;
    F77_CALL(dtrsm)("L", "U", "N", "N", &df_i, &df_i, &one_d, a_dd, &df_i,
        bsqrtinv_dd + bsqrtinv_idx_i, &df_i
        FCONE FCONE FCONE FCONE);

    bsqrtinv_idx_i += df_sq_i;
  }
  Free(corge_dd);
  Free(x_dd);
  Free(qraux_dd);
  Free(pivot_ii);
  Free(work_dd);
  Free(a_dd);

  /*
    foo <- function(ind) {
      bar <- crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind])
      baz <- matrix(BsqrtInv[,ind], nrow=df, ncol=df)
      crossprod(baz, bar)
    }
    qux <- lapply(1:V, foo)
  */
  statimg = PROTECT(allocMatrix(REALSXP, df_i, ncol_i));
  statimg_dd = REAL(statimg);

  x_dd = Calloc(df_i, double);
  for (col_i = 0; col_i < ncol_i; col_i++) {
    /* bar <- crossprod(sqrtSigma$X1res, sqrtSigma$res[,ind]) */
    one_d = 1.0;
    zero_d = 0.0;
    one_i = 1;

    F77_CALL(dgemv)("T", &n_i, &df_i, &one_d, x1res_dd, &n_i,
        res_dd + (col_i * n_i), &one_i, &zero_d, x_dd, &one_i FCONE); /* SNV replaced nrow_i with n_i */

    /*
      baz <- matrix(BsqrtInv[,ind], nrow=df, ncol=df)
      crossprod(baz, bar)
    */
    one_d = 1.0;
    zero_d = 0.0;
    one_i = 1;

    F77_CALL(dgemv)("T", &df_i, &df_i, &one_d, bsqrtinv_dd + (col_i * df_sq_i),
        &df_i, x_dd, &one_i, &zero_d, statimg_dd + (col_i * df_i), &one_i FCONE);
  }
  Free(x_dd);
  Free(bsqrtinv_dd);
  UNPROTECT(1); /* statimg */

  return statimg;
}

static const R_CallMethodDef callMethods[] = {
   {"pbj_pbjBootRobustX", (DL_FUNC) &pbj_pbjBootRobustX, 6},
   {NULL, NULL, 0}
};

void R_init_pbj(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_pbj(DllInfo *info) {
}
