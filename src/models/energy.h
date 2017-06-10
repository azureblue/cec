#ifndef ENERGY_H
#define	ENERGY_H
#include "../matrix.h"
#include "../errors.h"
#include "../cross_entropy_context.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#ifndef M_E
#define M_E 2.7182818284590452354
#endif

static double ZERO_EPSILON = 1.0e-32;

static inline double handle_zero(double d)
{
    if (d < ZERO_EPSILON)
        return ZERO_EPSILON;
    return d;
}

static inline double handle_cholesky_nan(double d)
{
    if (isnan(d))
        return handle_zero(0);
    return handle_zero(d);
}

/*
 * Cross-entropy function.
 */
typedef double (*cross_entropy_function) (struct cross_entropy_context *, const struct cec_matrix *);

/*
 * Implementations of cross-entropy function with respect to 
 * the Gaussian density family.
 */

#endif	/* ENERGY_H */
