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

SEXP pbj_sweep1Mult(SEXP res, SEXP stats) {
  SEXP resDim, statsDim;
  int *resDim_i, nrow, ncol, row, col, i;
  double *res_d, *stats_d;

  /* Type checking */
  if (isReal(res) == FALSE) {
    error("res must be a real 2d matrix");
  }

  resDim = getAttrib(res, R_DimSymbol);
  if (resDim == R_NilValue) {
    error("res must be a real 2d matrix");
  } else if (length(resDim) != 2) {
    error("res must be a real 2d matrix");
  }
  resDim_i = INTEGER(resDim);
  nrow = resDim_i[0];
  ncol = resDim_i[1];

  if (isReal(stats) == FALSE) {
    error("stats must be a real 1d vector");
  }

  statsDim = getAttrib(stats, R_DimSymbol);
  if (statsDim != R_NilValue) {
    error("stats must be a real 1d vector");
  }
  if (length(stats) != nrow) {
    error("stats must have the same length as the number of rows in res");
  }

  /* Reproduce this:
   *
   *  sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
   */

  /* NOTE: R stores matrices in column major order */
  res_d = REAL(res);
  stats_d = REAL(stats);
  i = 0;
  for (col = 0; col < ncol; col++ ) {
    for (row = 0; row < nrow; row++) {
      res_d[i] = res_d[i] * stats_d[row];
      i++;
    }
  }
  return res;
}

SEXP pbj_pbjBootRobustX(SEXP qr, SEXP res, SEXP x1res, SEXP h, SEXP df) {
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

  SEXP qr_qr, qr_qraux, elt, attr, dim, rsd, sweep;
  int *dim_ii, k_i, n_i, ny_i, res_nrow_i, res_ncol_i, row_i, col_i, rsd_idx_i,
      x1res_idx_i, sweep_idx_i, df_i, layer_i;
  double *rsd_dd, *h_dd, *x1res_dd, *sweep_dd;

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
      res_nrow_i = dim_ii[0];
      if (res_nrow_i != n_i) {
        error("qr and res must have the same number of rows");
      }

      /* ncol(res) */
      res_ncol_i = ny_i = dim_ii[1];
    } else {
      error("res must be a real 2d matrix");
    }
  } else {
    error("res must be a real 2d matrix");
  }

  /* Type checking for df */
  if (!isInteger(df) || length(df) != 1) {
    error("df must be an integer vector of length 1");
  }
  df_i = INTEGER(df)[0];

  /* Type checking for x1res */
  if (isReal(x1res)) {
    dim = getAttrib(x1res, R_DimSymbol);
    if (isInteger(dim) && length(dim) == 2) {
      dim_ii = INTEGER(dim);

      /* check nrow(x1res) == nrow(res) */
      if (dim_ii[0] != res_nrow_i) {
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

  /* Allocate space for rsd vector */
  rsd = PROTECT(allocVector(REALSXP, length(res)));
  dim = PROTECT(allocVector(INTSXP, 2));
  dim_ii = INTEGER(dim);
  dim_ii[0] = res_nrow_i;
  dim_ii[1] = res_ncol_i;
  setAttrib(rsd, R_DimSymbol, dim);
  UNPROTECT(1); /* dim */

  /* rsd <- qr.resid(sqrtSigma$QR, sqrtSigma$res) */
  rsd_dd = REAL(rsd);
  F77_NAME(dqrrsd)(REAL(qr_qr), &n_i, &k_i, REAL(qr_qraux), REAL(res), &ny_i, rsd_dd);

  /* rsd <- sweep(rsd, 1, 1-h, '/') */
  h_dd = REAL(h);
  rsd_idx_i = 0;
  for (col_i = 0; col_i < res_ncol_i; col_i++) {
    for (row_i = 0; row_i < res_nrow_i; row_i++) {
      rsd_dd[rsd_idx_i] = rsd_dd[rsd_idx_i] / (1 - h_dd[row_i]);
      rsd_idx_i++;
    }
  }

  /* Allocate space for sweep vector */
  sweep = PROTECT(allocVector(REALSXP, length(res) * df_i));
  dim = PROTECT(allocVector(INTSXP, 3));
  dim_ii = INTEGER(dim);
  dim_ii[0] = res_nrow_i;
  dim_ii[1] = res_ncol_i;
  dim_ii[2] = df_i;
  setAttrib(sweep, R_DimSymbol, dim);
  UNPROTECT(1); /* dim */

  /* baz <- rep(list(rsd), df)
     qux <- simplify2array(baz)
     corge <- sweep(qux, c(1,3), sqrtSigma$X1res, '*') */
  x1res_dd = REAL(x1res);
  x1res_idx_i = 0;
  sweep_dd = REAL(sweep);
  sweep_idx_i = 0;
  for (layer_i = 0; layer_i < df_i; layer_i++) {
    rsd_idx_i = 0;

    for (col_i = 0; col_i < res_ncol_i; ) {
      for (row_i = 0; row_i < res_nrow_i; row_i++) {
        sweep_dd[sweep_idx_i] = rsd_dd[rsd_idx_i] * x1res_dd[x1res_idx_i];
        sweep_idx_i++;
        rsd_idx_i++;
        x1res_idx_i++;
      }

      col_i++;

      /* Reset x1res index back to the first row in its current column, unless
       * we're about to break out of this inner loop */
      if (col_i < res_ncol_i) {
        x1res_idx_i -= res_nrow_i;
      }
    }
  }

  UNPROTECT(2); /* rsd, sweep */

  return sweep;
}

/* For debugging purposes in gdb */
SEXP pbj_inspect(SEXP x) {
  return x;
}

static const R_CallMethodDef callMethods[] = {
   {"pbj_pbjBootRobustX", (DL_FUNC) &pbj_pbjBootRobustX, 5},
   {"pbj_inspect", (DL_FUNC) &pbj_inspect, 1},
   {NULL, NULL, 0}
};

void R_init_pbj(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_pbj(DllInfo *info) {
}
