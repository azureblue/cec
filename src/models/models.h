#ifndef CEC_MODELS_H
#define CEC_MODELS_H

#include "../matrix.h"
#include "../cec_r.h"

struct cec_model * cec_create_ce_ctx_all(int n);
struct cec_model * cec_create_ce_ctx_covariance(const cec_mat * cov,
                                                const cec_mat * cov_inv);
struct cec_model * cec_create_ce_ctx_eigenvalues(int n, const double * given_evals);
struct cec_model * cec_create_ce_ctx_fixed_r(double r);
struct cec_model * cec_create_ce_ctx_spherical();
struct cec_model * cec_create_ce_ctx_diagonal();

struct cec_model * create_model(struct cec_model_spec * model_spec);

#endif //CEC_MODELS_H
