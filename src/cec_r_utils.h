#ifndef CEC_R_UTILS_H
#define	CEC_R_UTILS_H

#include <Rinternals.h>

#include "matrix.h"

/*
 * Utility functions, convertions between R matrix and cec_matrix structure.
 */
void copy_from_R_matrix(SEXP R_ma, struct cec_matrix * ma);

struct cec_matrix * create_from_R_matrix(SEXP R_matrix);

SEXP create_R_matrix(struct cec_matrix * m);

#endif	/* CEC_R_UTILS_H */

