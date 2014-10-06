
#ifndef MATRIX_UTILS_H
#define	MATRIX_UTILS_H

#include "matrix.h"

void cec_vector_outer_product(const double *restrict vec,
	struct cec_matrix *restrict output_matrix, int n);

double cec_matrix_trace(const struct cec_matrix * m);

void cec_matrix_multiply_sq(const struct cec_matrix *restrict m1,
	const struct cec_matrix *restrict m2, struct cec_matrix *restrict dest);

double cec_diagonal_product(const struct cec_matrix * m);

double cec_cov_cholesky_det(const struct cec_matrix * m,
	struct cec_matrix * temp);

int cec_cov_eigenvalues(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix, struct cec_matrix * workspace, double * values);

int cec_cov_cholesky(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix);

#endif	/* MATRIX_UTILS_H */

