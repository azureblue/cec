#include "energy_function.h"
#include "cov_utils.h"
#include <float.h>

static double _ZERO = DBL_MIN;
static double ZERO_DELTA = 1.11e-16;

static inline double handle_zero(double d) 
{
    if (d <= ZERO_DELTA) 
	d = _ZERO;
    return d;
}

static inline double det_z(double det) 
{
    if(isnan(det)) 
	return handle_zero(0);
    return handle_zero(det);
		
}

energy_function energy_function_for(enum density_family family)
{
    switch (family)
    {
	case ALL:
	    return h_all;
	case SPHERICAL:
	    return h_spherical;
	case DIAGONAL:
	    return h_diagonal;
	case FIXED_R:
	    return h_fixed_r;
	case GIVEN_COVARIANCE:
	    return h_given_covariance;
	case FIXEDEIGENVALUES:
	    return h_fixedeigenvalues;
    }

    return NULL;
}

double h_given_covariance(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{
    struct context_gc * cgc = (struct context_gc *) context->custom_context;
    struct cec_matrix * temp_matrix = context->temp_matrix;
    int n = context->n;

    cec_cov_multiply(cgc->i_given_cov, cov, temp_matrix);
    double trace = cec_cov_trace(temp_matrix);
    trace = handle_zero(trace);
    if (isnan(trace))
    {
	*context->last_error = POSITIVE_DEFINITE_ERROR;
	return NAN;
    }

    double det = cec_cov_cholesky_det(cgc->given_cov, temp_matrix);
    det = det_z(det);
    if (isnan(det))
    {
	*context->last_error = POSITIVE_DEFINITE_ERROR;
	return NAN;
    }

    return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * trace
	    + (1.0 / 2.0) * log(det);
}

double h_spherical(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{
    int n = context->n;

    double trace = cec_cov_trace(cov);
    trace = handle_zero(trace);
    if (isnan(trace))
    {
	*context->last_error = POSITIVE_DEFINITE_ERROR;
	return NAN;
    }
    
    return (n / 2.0) * log(2.0 * M_PI * M_E / n) + (n / 2.0) * log(trace);
}

double h_fixed_r(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{
    struct context_r * cr = (struct context_r *) context->custom_context;
    int n = context->n;
    double r = cr->r;

    return (n / 2.0) * log(2.0 * M_PI)
	    + (1.0 / (2.0 * r)) * cec_cov_trace(cov) + (n / 2.0) * log(r);
}

double h_diagonal(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{

    int n = context->n;

    double diagonal_product = cec_cov_diagonal_product(cov);
    diagonal_product = handle_zero(diagonal_product);
    if (isnan(diagonal_product))
    {
	*context->last_error = POSITIVE_DEFINITE_ERROR;
	return NAN;
    }

    return (n / 2.0) * log(2.0 * M_PI * M_E)
	    + (1.0 / 2.0) * log(diagonal_product);
}

double h_fixedeigenvalues(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{
    struct context_fe * cfe = (struct context_fe *) context->custom_context;
    struct cec_matrix * temp_matrix = context->temp_matrix;
    double * given_evals = cfe -> given_evals;
    double * evals	 = cfe -> evals;
    int n = context->n;

    int error = cec_cov_eigenvalues(cov, temp_matrix, cfe->workspace, evals);
    if (error != 0) 
    {
	*context->last_error = error;
	return NAN;
    }

    double e_sum = 0;
    for (int i = 0; i < n; i++)
	e_sum += evals[i] / given_evals[i];

    return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * e_sum + (1.0 / 2.0) * log(cfe->given_evals_product);

}

double h_all(const struct energy_function_context * context,
	const struct cec_matrix * cov)
{

    struct cec_matrix * temp_matrix = context->temp_matrix;
    int n = context->n;

    double det = cec_cov_cholesky_det(cov, temp_matrix);
    det = det_z(det);
    
    if (isnan(det))
    {
	*context->last_error = POSITIVE_DEFINITE_ERROR;
	return NAN;
    }
    return (n / 2.0) * log(2.0 * M_PI * M_E) + (1.0 / 2.0) * log(det);
}
