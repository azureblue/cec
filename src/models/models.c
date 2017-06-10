#include "models.h"

struct cec_model * create_model(struct cec_model_spec * model_spec)
{
    switch (model_spec->type)
    {
        case ALL:
            return cec_create_ce_ctx_all(model_spec->n);
        case SPHERICAL:
            return cec_create_ce_ctx_spherical();
        case DIAGONAL:
            return cec_create_ce_ctx_diagonal();
        case FIXED_R:
            return cec_create_ce_ctx_fixed_r((model_r_params(model_spec)->r));
        case GIVEN_COVARIANCE:
            return cec_create_ce_ctx_covariance(model_covariances_params(model_spec)->cov,
                                                model_covariances_params(model_spec)->cov_inv);
        case FIXEDEIGENVALUES:
            return cec_create_ce_ctx_eigenvalues(
                    model_eigenvalues_params(model_spec)->given_eigenvalues->len,
                    model_eigenvalues_params(model_spec)->given_eigenvalues->ar);
    }
    return NULL;
}
