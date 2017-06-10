#include "energy.h"
#include "../cov_utils.h"
#include "../model.h"

struct context_all {
    cec_mat * temp_matrix;
};

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov){

    struct cec_matrix * temp_matrix = ((struct context_all *)context->specific_context)->temp_matrix;
    int n = cov->n;

    double det = cec_cov_cholesky_det(cov, temp_matrix);
    det = handle_cholesky_nan(det);

    if (isnan(det)) {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }
    return (n / 2.0) * log(2.0 * M_PI * M_E) + (1.0 / 2.0) * log(det);
}

struct cec_model * cec_create_ce_ctx_all(int n) {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    struct context_all * c_all = alloc(struct context_all);
    c_all->temp_matrix = cec_matrix_create(n, n);
    model->cross_entropy_context->specific_context = c_all;
    model->cross_entropy = cross_entropy;
    return model;
}
