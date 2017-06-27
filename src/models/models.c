#include "models.h"

struct cec_model * create_model(struct cec_model_spec * model_spec)
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
            return cec_create_model_fixed_r(model_r_params(model_spec)->r);
        case COVARIANCE:
            return cec_create_model_covariance(model_covariances_params(model_spec)->cov,
                                               model_covariances_params(model_spec)->cov_inv);
        case EIGENVALUES:
            return cec_create_model_eigenvalues(
                    model_eigenvalues_params(model_spec)->given_eigenvalues->len,
                    model_eigenvalues_params(model_spec)->given_eigenvalues->ar);
    }
    return NULL;
}
