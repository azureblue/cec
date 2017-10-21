#include "models.h"

struct cec_model * create_model(cec_model_spec * model_spec)
{
    switch (model_spec->type)
    {
        case ALL:
            return cec_create_model_all(model_spec->n);
        case SPHERICAL:
            return cec_create_model_spherical();
        case DIAGONAL:
            return cec_create_model_diagonal();
        case FIXED_R:
            return cec_create_model_fixed_r(model_spec->r_params.r);
        case COVARIANCE:
            return cec_create_model_covariance(model_spec->covariances_params.cov,
                                               model_spec->covariances_params.cov_inv);
        case EIGENVALUES:
            return cec_create_model_eigenvalues(
                    model_spec->eigenvalues_params.given_eigenvalues->len,
                    model_spec->eigenvalues_params.given_eigenvalues->ar);
    }
    return NULL;
}
