#ifndef CENTERS_INIT_R_H
#define	CENTERS_INIT_R_H

#include <Rinternals.h>

/*
 * Center initialization methods - entry point functions called form R.
 */
SEXP init_kmeanspp_r(SEXP x, SEXP k);

SEXP init_random_r(SEXP x, SEXP k);

#endif	/* CENTERS_INIT_R_H */

