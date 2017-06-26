#include "cec_params_r.h"
#include "cec_r_utils.h"
#include "error_r.h"

cec_centers_par * get_centers_param(SEXP centers_param_r) {
    cec_centers_par *centers = alloc(cec_centers_par);
    if (!parse_init_method(CHAR(STRING_ELT(get_named_element(centers_param_r, "init.method"), 0)), &centers->init_m))
        error_r(LIBRARY_DEFECT_ERROR);
    SEXP var_centers = get_named_element(centers_param_r, "var.centers");
    centers->var_centers = vec_i_create_from(LENGTH(var_centers), INTEGER(var_centers));
    centers->centers_mat = NULL;
    if (centers->init_m == NONE)
        centers->centers_mat = create_from_R_matrix(get_named_element(centers_param_r, "mat"));
    return centers;
}

cec_control_par * get_control_param(SEXP control_param_r) {
    cec_control_par *control = alloc(cec_control_par);
    control->min_card = asInteger(get_named_element(control_param_r, "min.card"));
    control->max_iterations = asInteger(get_named_element(control_param_r, "max.iters"));
    control->starts = asInteger(get_named_element(control_param_r, "starts"));
    control->threads = asInteger(get_named_element(control_param_r, "threads"));
    return control;
}

cec_models_par * get_models_param(SEXP models_param_r, int n) {
    int len = LENGTH(models_param_r);
    cec_models_par *model_par = alloc_fam(cec_models_par, struct cec_model_spec, len);
    model_par->len = len;
    for (int i = 0; i < len; i++) {
        SEXP model_r = VECTOR_ELT(models_param_r, i);
        enum density_family type = (enum density_family) asInteger(get_named_element(model_r, "type"));
        SEXP params = get_named_element(model_r, "params");
        struct cec_model_spec *spec = &model_par->model_specs[i];
        spec->n = n;
        spec->type = type;
        switch (type) {
            case FIXED_R:
                spec->type_specific_params = alloc(struct cec_model_r_params);
                ((struct cec_model_r_params*) spec->type_specific_params)->r = asReal(get_named_element(params, "r"));
                break;
            case GIVEN_COVARIANCE: {
                struct cec_model_covariances_params *cov_params = spec->type_specific_params = alloc(
                struct cec_model_covariances_params);
                cov_params->cov = create_from_R_matrix(get_named_element(params, "cov"));
                cov_params->cov_inv = create_from_R_matrix(get_named_element(params, "cov.inv"));
                break;
            }
            case FIXEDEIGENVALUES: {
                struct cec_model_eigenvalues_params *eigen_params = spec->type_specific_params = alloc(
                struct cec_model_eigenvalues_params);
                SEXP eigenvalues = get_named_element(params, "eigenvalues");
                eigen_params->given_eigenvalues = vec_d_create_from(LENGTH(eigenvalues), REAL(eigenvalues));
                break;
            }
            case SPHERICAL:break;
            case DIAGONAL:break;
            case ALL:break;
        }
    }
    return model_par;
}
