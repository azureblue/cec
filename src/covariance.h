#ifndef COVARIANCE_MATRIX_H
#define	COVARIANCE_MATRIX_H

#include "matrix.h"

void
covariance_add_point(const struct cec_matrix * covariance,
		struct cec_matrix * new_covarioance, const double * mean,
		double const * point, int card, struct cec_matrix * t_matrix);

void
covariance_remove_point(const struct cec_matrix * covariance,
		struct cec_matrix * new_covarioance, const double * mean,
		double const * point, int card, struct cec_matrix * t_matrix);

#endif	/* COVARIANCE_MATRIX_H */

