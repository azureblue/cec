#ifndef CEC_R_H
#define	CEC_R_H

#include <Rinternals.h>

/*
 * Entry point to CEC - called from R.
 */
SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r);

#endif	/* CEC_R_H */
