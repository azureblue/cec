
#ifndef MATRIX_UTILS_H
#define	MATRIX_UTILS_H

#include "matrix.h"

int cec_eigenvalues_sm(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix, double * work, int lwork, double * values);

#endif	/* MATRIX_UTILS_H */

