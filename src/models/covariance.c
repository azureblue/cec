#include "energy.h"
#include "../cov_utils.h"
#include "../alloc.h"
#include "../model.h"

struct context_gc
{
    struct cec_matrix * given_cov;
    struct cec_matrix * i_given_cov;
    struct cec_matrix * temp_matrix;
};

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov) {
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

struct cec_model * cec_create_ce_ctx_covariance(const struct cec_matrix * cov,
                                         const struct cec_matrix * cov_inv) {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    struct context_gc * c_gc =  alloc(struct context_gc);
    int n = cov->n;
    c_gc->given_cov = cec_matrix_create_copy(cov);
    c_gc->i_given_cov = cec_matrix_create_copy(cov_inv);
    c_gc->temp_matrix = cec_matrix_create(n, n);
    model->cross_entropy_context->specific_context = c_gc;
    model->cross_entropy = cross_entropy;
    return model;
}

