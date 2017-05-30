#include <R.h>
#include <Rdefines.h>
#include "alloc.h"
#include "cec_r.h"
#include "cec_r_utils.h"
#include "cec.h"
#include "error_r.h"
#include "alloc_utils_r.h"

static struct cec_model ** create_cec_models(SEXP type, SEXP params, int m, int k, int n);
static struct cec_model * create_model_from_R_params(enum density_family family, SEXP param, int n);
static SEXP create_R_result(cec_out* , int m);

static cec_centers * get_centers_param(SEXP centers_param_r) {
    cec_centers *centers = alloc(cec_centers);
    if (!parse_init_method(CHAR(STRING_ELT(get_named_element(centers_param_r, "init.method"), 0)), &centers->init_m))
        error_r(LIBRARY_DEFECT_ERROR);
    SEXP var_centers = get_named_element(centers_param_r, "var.centers");
    centers->var_centers = vec_i_create_from(LENGTH(var_centers), INTEGER(var_centers));
    centers->centers_mat = NULL;
    if (centers->init_m == NONE)
        centers->centers_mat = create_from_R_matrix(get_named_element(centers_param_r, "mat"));
    return centers;
}

static cec_control * get_control_param(SEXP control_param_r) {
    cec_control *control = alloc(cec_control);
    control->min_card = asInteger(get_named_element(control_param_r, "min.card"));
    control->max_iterations = asInteger(get_named_element(control_param_r, "max.iters"));
    control->starts = asInteger(get_named_element(control_param_r, "starts"));
    return control;
}

SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP type, SEXP params)
{
    init_mem_mg(error_r_mem_error);
    int m = Rf_nrows(x);
    int n = Rf_ncols(x);
    cec_centers *centers = get_centers_param(centers_param_r);
    cec_control *control = get_control_param(control_param_r);
    int vc_len = centers->var_centers->len;
    int vc_max_k = -1;
    for (int i = 0; i < vc_len; i++)
        vc_max_k = vc_max_k < centers->var_centers->ar[i] ? centers->var_centers->ar[i] : vc_max_k;
    cec_mat * x_mat = create_from_R_matrix(x);
    cec_mat * c_mat = cec_matrix_create(vc_max_k, n);
    if (centers->centers_mat)
        cec_matrix_copy_data(centers->centers_mat, c_mat);

    struct cec_model ** cec_models = create_cec_models(type, params, m, vc_max_k, n);

    double best_energy = BIG_DOUBLE;
    cec_out *best_result = create_cec_output(m, vc_max_k, n, control->max_iterations);
    res_code all_res = UNKNOWN_ERROR;

    int starts = control->starts;

    for (int start = 0; start < starts; start++) {
        for (int vc = 0; vc < vc_len; vc++) {
            int vc_k = centers->var_centers->ar[vc];
            m_state ms = m_current_state();
            c_mat->m = vc_k;
            cec_init_centers(x_mat, c_mat, centers->init_m);
            cec_in *in = create_cec_input(x_mat, c_mat, cec_models, control->max_iterations, control->min_card);
            cec_out *out = create_cec_output(m, vc_k, n, control->max_iterations);
            cec_ctx *ctx = create_cec_context(in, out);
            int res = cec_start(ctx);
            if (res == NO_ERROR) {
                all_res = NO_ERROR;
                double energy = cec_final_energy(out);
                if (energy < best_energy) {
                    cec_copy_results_content(out, best_result);
                    best_energy = energy;
                }
            }
            m_reset_state(ms);
        }
    }

    SEXP result = NULL;

    if (all_res == NO_ERROR) {
        PROTECT(result = create_R_result(best_result, m));
        UNPROTECT(1);
        free_mem_mg();
        return result;
    }
    free_mem_mg();
    error_r(all_res);

    return NULL;
}

static SEXP create_R_result(cec_out * results, int m)
{
    int k = results->centers->m;
    int trimmed_size = results->iterations + 1;
    SEXP energy_vector;
    SEXP clusters_number_vector;
    SEXP assignment_vector;
    SEXP covariance_list;
    SEXP centers_matrix;
    SEXP iterations;

    vec_i *cluster_number_trimmed = vec_i_create_from(trimmed_size, results->clusters_number->ar);
    vec_d *energy_trimmed = vec_d_create_from(trimmed_size, results->energy->ar);

    PROTECT(energy_vector = allocVector(REALSXP, trimmed_size));
    PROTECT(clusters_number_vector = allocVector(INTSXP, trimmed_size));
    PROTECT(assignment_vector = allocVector(INTSXP, m));
    PROTECT(covariance_list = allocVector(VECSXP, k));
    PROTECT(iterations = allocVector(INTSXP, 1));
    PROTECT(centers_matrix = create_R_matrix(results->centers));

    INTEGER(iterations)[0] = results->iterations;


    vec_d_copy_to(energy_trimmed, REAL(energy_vector));
    vec_i_copy_to(cluster_number_trimmed, INTEGER(clusters_number_vector));
    vec_i_copy_to(results->clustering_vector, INTEGER(assignment_vector));
    for (int i = 0; i < m; i++)
        INTEGER(assignment_vector)[i]++;


    for (int i = 0; i < k; i++)
    {
        SEXP covariance;
        PROTECT(covariance = create_R_matrix(results->covriances->mats[i]));
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

static struct cec_model ** create_cec_models(SEXP type, SEXP params, int m, int k, int n)
{
    struct cec_model ** cec_models = alloc_n(struct cec_model *, k);
    for (int i = 0; i < k; i++)
    {
        SEXP param = VECTOR_ELT(params, i);
        enum density_family family = (enum density_family) INTEGER(type)[i];
        cec_models[i] = create_model_from_R_params(family, param, n);
    }
    return cec_models;
}

struct cec_model * create_model_from_R_params(enum density_family family, SEXP param, int n)
{
    struct cec_model *model = alloc(struct cec_model);
    switch (family)
    {
        case ALL:
            model->cross_entropy = h_all;
            model->cross_entropy_context = create_cross_entropy_context_all(n);
            break;
        case SPHERICAL:
            model->cross_entropy = h_spherical;
            model->cross_entropy_context = create_cross_entropy_context_spherical();
            break;
        case DIAGONAL:
            model->cross_entropy = h_diagonal;
            model->cross_entropy_context = create_cross_entropy_context_diagonal();
            break;
        case FIXED_R:
            model->cross_entropy = h_fixed_r;
            model->cross_entropy_context = create_cross_entropy_context_fixedr(asReal(param));
            break;
        case GIVEN_COVARIANCE: {
            struct cec_matrix *cov = create_from_R_matrix(VECTOR_ELT(param, 0));
            struct cec_matrix *cov_i = create_from_R_matrix(VECTOR_ELT(param, 1));
            model->cross_entropy = h_given_covariance;
            model->cross_entropy_context = create_cross_entropy_context_covariance(n, cov, cov_i);

            break;
        }
        case FIXEDEIGENVALUES:
            model->cross_entropy = h_fixedeigenvalues;
            model->cross_entropy_context = create_cross_entropy_context_eigenvalues(n, REAL(param));
    }

    return model;
}
