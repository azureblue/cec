#ifndef CEC_R_H
#define	CEC_R_H

#include <Rinternals.h>

/*
 * Entry point to CEC - called from R.
 */
SEXP cec_r(SEXP x, SEXP centers, SEXP iter_max, SEXP type, SEXP card_min,
	SEXP params);

#endif	/* CEC_R_H */

