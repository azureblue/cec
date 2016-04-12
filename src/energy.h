#ifndef ENERGY_H
#define	ENERGY_H

#include "matrix.h"
#include "errors.h"
#include "cross_entropy_context.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#ifndef M_E
#define M_E 2.7182818284590452354
#endif

/*
 * Cross-entropy function.
 */
typedef double (*cross_entropy_function) (struct cross_entropy_context *,
        const struct cec_matrix *);

/*
 * Implementations of cross-entropy function with respect to 
 * the Gaussian density family.
 */
double h_given_covariance
(struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_spherical
(struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_all
(struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_diagonal
(struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_fixed_r
(struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_fixedeigenvalues
(struct cross_entropy_context * context, const struct cec_matrix * cov);

#endif	/* ENERGY_H */
