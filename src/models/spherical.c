#include "energy.h"
#include "../cov_utils.h"
#include "../alloc.h"
#include "../model.h"

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov) {
    int n = cov->n;
    double trace = cec_cov_trace(cov);
    trace = handle_zero(trace);
    if (isnan(trace)) {
        context->last_error = INVALID_COVARIANCE_ERROR;
        return NAN;
    }
    return (n / 2.0) * log(2.0 * M_PI * M_E / n) + (n / 2.0) * log(trace);
}

struct cec_model * cec_create_model_spherical() {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    model->cross_entropy = cross_entropy;
    return model;
}

