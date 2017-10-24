#ifndef CEC_MODELS_H
#define CEC_MODELS_H

#include "../matrix.h"
#include "../cec_r.h"
#include "../cec_params.h"
#include "../model.h"

cec_model * cec_create_model_all(int n);
cec_model * cec_create_model_covariance(const cec_mat *cov,
                                               const cec_mat *cov_inv);
cec_model * cec_create_model_eigenvalues(int n, const double *given_evals);
cec_model * cec_create_model_fixed_r(double r);
cec_model * cec_create_model_spherical();
cec_model * cec_create_model_diagonal();

cec_model * create_model(cec_model_spec * model_spec);

#endif //CEC_MODELS_H
