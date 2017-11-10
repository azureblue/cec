#include "model.h"

static inline double energy(double cross_entropy, int card, int m)
{
    double p = (card / (double) m);
    return p * (-log(p) + cross_entropy);
}

double cluster_energy(const cec_model * model, cec_mat * cov, int card, int m)
{
    double ce = model->cross_entropy(model->cross_entropy_context, cov);
    return isnan(ce) ? ce : energy(ce, card, m);
}

res_code cluster_energy_get_last_error(const cec_model * model)
{
    return model->cross_entropy_context->last_error;
}
