#include "energy.h"
#include "../cov_utils.h"
#include "../alloc.h"
#include "../model.h"

struct context_r
{
    double r;
};

static double cross_entropy(struct cross_entropy_context * context, const struct cec_matrix * cov){
    struct context_r * cr = (struct context_r *) context->specific_context;
    int n = cov->n;
    double r = cr->r;

    return (n / 2.0) * log(2.0 * M_PI)
           + (1.0 / (2.0 * r)) * cec_cov_trace(cov) + (n / 2.0) * log(r);
}

struct cec_model * cec_create_ce_ctx_fixed_r(double r) {
    struct cec_model *model = alloc(struct cec_model);
    model->cross_entropy_context = alloc(cross_entropy_ctx);
    struct context_r * c_r =  alloc(struct context_r);
    c_r->r = r;
    model->cross_entropy_context->specific_context = c_r;
    model->cross_entropy = cross_entropy;
    return model;
}
