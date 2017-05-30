#include <R.h>
#include <Rdefines.h>
#include "alloc.h"
#include "cec_r.h"
#include "cec_r_utils.h"
#include "cec.h"
#include "error_r.h"
#include "init_utils_r.h"

static struct cec_model ** create_cec_models(SEXP type, SEXP params, int m, int k, int n);
static struct cec_model * create_model_from_R_params(enum density_family family, SEXP param, int n);
static SEXP create_R_result(cec_out*);
res_code cec_perform(cec_mat *x_mat, cec_centers_par *centers, cec_control_par *control, cec_models_par *cec_models,
                     cec_out **results);

static cec_centers_par * get_centers_param(SEXP centers_param_r) {
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

static cec_control_par * get_control_param(SEXP control_param_r) {
    cec_control_par *control = alloc(cec_control_par);
    control->min_card = asInteger(get_named_element(control_param_r, "min.card"));
    control->max_iterations = asInteger(get_named_element(control_param_r, "max.iters"));
    control->starts = asInteger(get_named_element(control_param_r, "starts"));
    return control;
}

static cec_models_par * get_models_param(SEXP models_param_r, int n) {
    int len = LENGTH(models_param_r);
    cec_models_par *models = alloc_fam(cec_models_par, struct cec_model_spec, len);
    models->len = len;
    for (int i = 0; i < len; i++) {
        SEXP model_r = VECTOR_ELT(models_param_r, i);
        enum density_family type = (enum density_family) asInteger(get_named_element(model_r, "type"));
        SEXP params = get_named_element(model_r, "params");
        struct cec_model_spec *spec = &models->model_specs[i];
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

    return models;
}

SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP models_param_r)
{
    cec_init_env();
    int m = Rf_nrows(x);
    int n = Rf_ncols(x);
    cec_centers_par *centers = get_centers_param(centers_param_r);
    cec_control_par *control = get_control_param(control_param_r);
    cec_models_par *models = get_models_param(models_param_r, n);
    int vc_len = centers->var_centers->len;
    int vc_max_k = -1;
    for (int i = 0; i < vc_len; i++)
        vc_max_k = vc_max_k < centers->var_centers->ar[i] ? centers->var_centers->ar[i] : vc_max_k;

    cec_out * out;
    res_code all_res = cec_perform(create_from_R_matrix(x), centers, control, models, &out);

    SEXP result = NULL;

    if (all_res == NO_ERROR) {
        PROTECT(result = create_R_result(out));
        UNPROTECT(1);
        cec_clean_env();
        return result;
    }
    cec_clean_env();
    error_r(all_res);

    return NULL;
}

static SEXP create_R_result(cec_out * out)
{
    int m = out->clustering_vector->len;
    int k = out->initial_k;
    int n = out->centers->n;
    int trimmed_size = out->iterations + 1;
    SEXP energy_vector;
    SEXP clusters_number_vector;
    SEXP assignment_vector;
    SEXP covariance_list;
    SEXP centers_matrix;
    SEXP iterations;

    vec_i *cluster_number_trimmed = vec_i_create_from(trimmed_size, out->clusters_number->ar);
    vec_d *energy_trimmed = vec_d_create_from(trimmed_size, out->energy->ar);
    cec_mat *trimmed_centers = cec_matrix_create(k, n);
    for (int i = 0; i < k; i++)
        array_copy(cec_matrix_const_row(out->centers, i), cec_matrix_row(trimmed_centers, i), n);

    PROTECT(energy_vector = allocVector(REALSXP, trimmed_size));
    PROTECT(clusters_number_vector = allocVector(INTSXP, trimmed_size));
    PROTECT(assignment_vector = allocVector(INTSXP, m));
    PROTECT(covariance_list = allocVector(VECSXP, k));
    PROTECT(iterations = allocVector(INTSXP, 1));
    PROTECT(centers_matrix = create_R_matrix(trimmed_centers));

    INTEGER(iterations)[0] = out->iterations;

    vec_d_copy_to(energy_trimmed, REAL(energy_vector));
    vec_i_copy_to(cluster_number_trimmed, INTEGER(clusters_number_vector));
    vec_i_copy_to(out->clustering_vector, INTEGER(assignment_vector));
    for (int i = 0; i < m; i++)
        INTEGER(assignment_vector)[i]++;


    for (int i = 0; i < k; i++)
    {
        SEXP covariance;
        PROTECT(covariance = create_R_matrix(out->covriances->mats[i]));
        SET_VECTOR_ELT(covariance_list, i, covariance);
    }

    SEXP ret;
    PROTECT(ret = allocList(6));
    
    SEXP ret_s = ret;
    
    SETCAR(ret, assignment_vector);
    SET_TAG(ret, install("cluster"));
    ret = CDR(ret);
    SETCAR(ret, centers_matrix);
    SET_TAG(ret, install("centers"));
    ret = CDR(ret);
    SETCAR(ret, energy_vector);
    SET_TAG(ret, install("energy"));
    ret = CDR(ret);
    SETCAR(ret, clusters_number_vector);
    SET_TAG(ret, install("nclusters"));
    ret = CDR(ret);
    SETCAR(ret, covariance_list);
    SET_TAG(ret, install("covariances"));
    ret = CDR(ret);
    SETCAR(ret, iterations);
    SET_TAG(ret, install("iterations"));
    
    UNPROTECT(k + 7);
    
    return ret_s;
}
