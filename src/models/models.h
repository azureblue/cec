#ifndef CEC_MODELS_H
#define CEC_MODELS_H

#include "../matrix.h"
#include "../cec_r.h"
#include "../cec_params.h"

struct cec_model * cec_create_model_all(int n);
struct cec_model * cec_create_model_covariance(const cec_mat *cov,
                                               const cec_mat *cov_inv);
struct cec_model * cec_create_model_eigenvalues(int n, const double *given_evals);
struct cec_model * cec_create_model_fixed_r(double r);
struct cec_model * cec_create_model_spherical();
struct cec_model * cec_create_model_diagonal();

struct cec_model * create_model(struct cec_model_specification * model_spec);

#endif //CEC_MODELS_H
