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
static SEXP create_R_result(cec_res* , int m);
static cec_centers * get_centers_param(SEXP centers_param_r, int n);
static cec_control * get_control_param(SEXP control_param_r);

SEXP cec_r(SEXP x, SEXP centers_param_r, SEXP control_param_r, SEXP type, SEXP params)
{
    init_mem_mg(error_r_mem_error);
    int m = Rf_nrows(x);
    int n = Rf_ncols(x);    
    cec_centers *centers = get_centers_param(centers_param_r, n);
    cec_control *control = get_control_param(control_param_r);
    cec_mat * c_mat = centers->centers_mat;
    int max_k = c_mat->m;

    cec_mat * x_mat = create_from_R_matrix(x);
    struct cec_model ** cec_models = create_cec_models(type, params, m, max_k, n);
    
    double best_energy = BIG_DOUBLE;
    cec_res *best_result = create_cec_result(m, max_k, n, control->max_iterations);
    int all_res = -1;
    
    int starts = control->starts;
    int vc_len = centers->var_centers_len;
    
    for (int start = 0; start < starts; start++) {
        for (int vc = 0; vc < vc_len; vc++) {
            int k = centers->var_centers[vc];
            m_state ms = m_current_state();
            c_mat->m = k;
            cec_init_centers(x_mat, c_mat, centers->init_m);
            cec_in *c_input = create_cec_input(x_mat, c_mat, cec_models, control->max_iterations, control->min_card);
            cec_res *c_res = create_cec_result(m, k, n, control->max_iterations);
            struct cec_context *cec_ctx = create_cec_context(c_input, c_res);
            int res = cec(cec_ctx);
            if (res == NO_ERROR) {
                all_res = NO_ERROR;
                double energy = cec_final_energy(c_res);
                if (energy < best_energy) {
                    cec_copy_results_content(c_res, best_result, m);
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

static SEXP create_R_result(cec_res * results, int m)
{
    int iters = results->iterations;
    int k = results->centers->m;
    int output_size = iters + 1;
    SEXP energy_vector;
    SEXP clusters_number_vector;
    SEXP assignment_vector;
    SEXP covariance_list;
    SEXP centers_matrix;
    SEXP iterations;

    PROTECT(energy_vector = allocVector(REALSXP, output_size));
    PROTECT(clusters_number_vector = allocVector(INTSXP, output_size));
    PROTECT(assignment_vector = allocVector(INTSXP, m));
    PROTECT(covariance_list = allocVector(VECSXP, k));
    PROTECT(iterations = allocVector(INTSXP, 1));
    PROTECT(centers_matrix = create_R_matrix(results->centers));

    INTEGER(iterations)[0] = iters;

    for (int i = 0; i < output_size; i++)
    {
        REAL(energy_vector)[i] = results->energy[i];
        INTEGER(clusters_number_vector)[i] = results->clusters_number[i];
    }

    for (int i = 0; i < m; i++)
        INTEGER(assignment_vector)[i] = results->clustering_vector[i] + 1;

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
        enum density_family family = INTEGER(type)[i];
        cec_models[i] = create_model_from_R_params(family, param, n);
    }
    return cec_models;
}

struct cec_model * create_model_from_R_params(enum density_family family, SEXP param, int n)
{
    struct cec_model * model = alloc(struct cec_model);
    model->family = family;
    
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
            case GIVEN_COVARIANCE:
            {
                struct cec_matrix * cov = create_from_R_matrix(VECTOR_ELT(param, 0));
                struct cec_matrix * cov_i = create_from_R_matrix(VECTOR_ELT(param, 1));
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

cec_centers * get_centers_param(SEXP centers_param_r, int n) {
    cec_centers *centers = alloc(cec_centers);
    if (!parse_init_method(CHAR(STRING_ELT(get_named_element(centers_param_r, "init.method"), 0)), &centers->init_m))
        error_r(LIBRARY_DEFECT_ERROR);
    SEXP var_centers = get_named_element(centers_param_r, "var.centers");
    centers->var_centers_len = LENGTH(var_centers);
    centers->var_centers = INTEGER(var_centers);
    int max_centers_num = centers->var_centers[0];
    for (int i = 1; i < centers->var_centers_len; i++)
        max_centers_num = max_centers_num < centers->var_centers[i] ? centers->var_centers[i] : max_centers_num;
    centers->centers_mat = cec_matrix_create(max_centers_num, n);
    if (centers->init_m == NONE)
        copy_from_R_matrix(get_named_element(centers_param_r, "mat"), centers->centers_mat);
    return centers;
}

static cec_control * get_control_param(SEXP control_param_r) {
    cec_control *control = alloc(cec_control);
    control->min_card = asInteger(get_named_element(control_param_r, "min.card"));
    control->max_iterations = asInteger(get_named_element(control_param_r, "max.iters"));
    control->starts = asInteger(get_named_element(control_param_r, "starts"));
    return control;
}
