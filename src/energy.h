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
 * Overall cluster energy (cost function).
 */
static inline double cluster_energy(const int m, const double cross_entropy,
	const int card)
{
    double p = (card / (double) m);
    return p * (-log(p) + cross_entropy);
}

/*
 * Internal energy function (cross-entropy).
 */
typedef double (*cross_entropy_function) (const struct cross_entropy_context *,
	const struct cec_matrix *);

cross_entropy_function cross_entropy_for(enum density_family family);

/*
 * Implementations of cross-entropy function with respect to 
 * the Gaussian density family.
 */
double h_given_covariance
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_spherical
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_all
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_diagonal
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_fixed_r
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

double h_fixedeigenvalues
(const struct cross_entropy_context * context, const struct cec_matrix * cov);

#endif	/* ENERGY_H */
