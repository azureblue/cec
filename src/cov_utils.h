#ifndef COV_UTILS_H
#define	COV_UTILS_H

#include <stdbool.h>
#include "matrix.h"

/*
 * Set of functions operating on covariance matrix.
 */

double cec_cov_trace(const struct cec_matrix * m);

double cec_cov_diagonal_product(const struct cec_matrix * m);

void cec_cov_multiply(const struct cec_matrix *restrict m1,
        const struct cec_matrix *restrict m2, struct cec_matrix *restrict dest);

double cec_cov_cholesky_det(const struct cec_matrix * m,
        struct cec_matrix * temp);

bool cec_cov_eigenvalues(const struct cec_matrix * sym_matrix,
                         struct cec_matrix * temp_matrix, struct cec_matrix * workspace, double * values);

bool cec_cov_cholesky(const struct cec_matrix * sym_matrix,
                      struct cec_matrix * temp_matrix);

void cec_cov_add_point(struct cec_matrix *restrict dest_covariance, const struct cec_matrix *restrict covariance, 
        const double *restrict mean, double const *restrict point, int card);

void cec_cov_remove_point(struct cec_matrix *restrict dest_covariance, const struct cec_matrix *restrict covariance, 
        const double *restrict mean, double const *restrict point, int card);

#endif	/* COV_UTILS_H */
