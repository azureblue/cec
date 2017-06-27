#include "energy.h"
#include "../cov_utils.h"
#include "../alloc.h"
#include "../model.h"

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov) {
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

struct cec_model * cec_create_model_diagonal() {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    model->cross_entropy = cross_entropy;
    return model;
}
