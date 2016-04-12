#include "energy.h"
#include "cov_utils.h"
#include <float.h>

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

double h_given_covariance(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{
    struct context_gc * cgc = (struct context_gc *) context->specific_context;
    struct cec_matrix * temp_matrix = cgc->temp_matrix;
    int n = cov->n;

    cec_cov_multiply(cgc->i_given_cov, cov, temp_matrix);
    double trace = cec_cov_trace(temp_matrix);
    trace = handle_zero(trace);
    double det = cec_cov_cholesky_det(cgc->given_cov, temp_matrix);
    det = handle_cholesky_nan(det);
    if (isnan(trace) || isnan(det))
    {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }

    return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * trace
            + (1.0 / 2.0) * log(det);
}

double h_spherical(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{
    int n = cov->n;

    double trace = cec_cov_trace(cov);
    trace = handle_zero(trace);
    if (isnan(trace))
    {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }

    return (n / 2.0) * log(2.0 * M_PI * M_E / n) + (n / 2.0) * log(trace);
}

double h_fixed_r(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{
    struct context_r * cr = (struct context_r *) context->specific_context;
    int n = cov->n;
    double r = cr->r;

    return (n / 2.0) * log(2.0 * M_PI)
            + (1.0 / (2.0 * r)) * cec_cov_trace(cov) + (n / 2.0) * log(r);
}

double h_diagonal(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{

    int n = cov->n;

    double diagonal_product = cec_cov_diagonal_product(cov);
    diagonal_product = handle_zero(diagonal_product);
    if (isnan(diagonal_product))
    {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }

    return (n / 2.0) * log(2.0 * M_PI * M_E)
            + (1.0 / 2.0) * log(diagonal_product);
}

double h_fixedeigenvalues(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{
    struct context_fe * cfe = (struct context_fe *) context->specific_context;
    struct cec_matrix * temp_matrix = cfe->temp_matrix;
    double * given_evals = cfe -> evals_given;
    double * evals = cfe -> evals_temp;
    int n = cov->n;

    int error = cec_cov_eigenvalues(cov, temp_matrix, cfe->workspace, evals);
    if (error)
    {
        context->last_error = error;
        return NAN;
    }

    double e_sum = 0;
    for (int i = 0; i < n; i++)
        e_sum += evals[i] / given_evals[i];

    return (n / 2.0) * log(2.0 * M_PI) + (1.0 / 2.0) * e_sum + (1.0 / 2.0) * log(cfe->given_evals_product);
}

double h_all(struct cross_entropy_context * context,
        const struct cec_matrix * cov)
{

    struct cec_matrix * temp_matrix = ((struct context_all *)context->specific_context)->temp_matrix;
    int n = cov->n;

    double det = cec_cov_cholesky_det(cov, temp_matrix);
    det = handle_cholesky_nan(det);

    if (isnan(det))
    {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }
    return (n / 2.0) * log(2.0 * M_PI * M_E) + (1.0 / 2.0) * log(det);
}
