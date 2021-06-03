#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP pbjBootRobustT(SEXP res, SEXP stats) {
  SEXP resDim, statsDim;
  int *resDim_i, nrow, ncol, row, col, i;
  double *res_d, *stats_d;

  /* Type checking */
  if (isReal(res) == FALSE) {
    error("res must be numeric");
  }

  resDim = getAttrib(res, R_DimSymbol);
  if (resDim == R_NilValue) {
    error("res must be a numeric 2d matrix");
  } else if (length(resDim) != 2) {
    error("res must be a numeric 2d matrix");
  }
  resDim_i = INTEGER(resDim);
  nrow = resDim_i[0];
  ncol = resDim_i[1];

  if (isReal(stats) == FALSE) {
    error("stats must be numeric");
  }

  statsDim = getAttrib(stats, R_DimSymbol);
  if (statsDim != R_NilValue) {
    error("stats must be a numeric 1d vector");
  }
  if (length(stats) != nrow) {
    error("stats must have the same length as the number of rows in res");
  }

  /* Reproduce this:
   *   sqrtSigma$res = sweep(sqrtSigma$res, 1, rboot(n)/sqrt(1-h), '*')
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

static const R_CallMethodDef callMethods[] = {
   {"pbj_pbjBootRobustT", (DL_FUNC) &pbjBootRobustT, 2},
   {NULL, NULL, 0}
};

void R_init_pbj(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_pbj(DllInfo *info) {
}
