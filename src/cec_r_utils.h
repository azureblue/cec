#ifndef CECR_UTILS_H
#define	CECR_UTILS_H

#include <Rinternals.h>

#include "matrix.h"

void copy_from_R_matrix(SEXP R_ma, struct cec_matrix * ma);

struct cec_matrix * create_from_R_matrix(SEXP R_matrix);

SEXP create_R_matrix(struct cec_matrix * m);

#endif	/* CECR_UTILS_H */

