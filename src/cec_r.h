#ifndef CECR_H
#define	CECR_H

#include <Rinternals.h>

SEXP cec_r(SEXP x, SEXP centers, SEXP iter_max, SEXP type, SEXP card_min,
	SEXP params);

#endif	/* CECR_H */

