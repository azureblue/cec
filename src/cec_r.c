#include <R.h>
#include <Rdefines.h>
#include "alloc.h"
#include "cec_r.h"
#include "cec_r_utils.h"
#include "cec.h"
#include "error_r.h"
#include "model.h"
#include "rand.h"
#include "alloc_utils_r.h"

static struct cec_model ** create_cec_models(SEXP type, SEXP params, int m, int k, int n);

static struct cec_model * create_model_from_R_params(enum density_family family, SEXP param, int n);

static SEXP create_R_result(struct cec_context * cec_c, int m, int k);

SEXP cec_r(SEXP x, SEXP centers, SEXP iter_max, SEXP type, SEXP card_min, SEXP params)
{
    init_mem_mg(error_r_mem_error);

    int m = Rf_nrows(x);
    int k = Rf_nrows(centers);
    int n = Rf_ncols(x);
    int iteration_max = Rf_asInteger(iter_max);
    int card_min_int = Rf_asInteger(card_min);

    struct cec_matrix * X = create_from_R_matrix(x);
    struct cec_matrix * C = create_from_R_matrix(centers);
    struct cec_model ** cec_models = create_cec_models(type, params, m, k, n);
    struct cec_context * cec_ctx = create_cec_context(X, C, cec_models, iteration_max, card_min_int);
    int res = cec(cec_ctx);

    SEXP result = NULL;

    if (res == NO_ERROR) {
        PROTECT(result = create_R_result(cec_ctx, m, k));
        UNPROTECT(1);
        free_mem_mg();
        return result;
    }
    free_mem_mg();
    error_r(res);

    return NULL;
}

static SEXP create_R_result(struct cec_context * cec_c, int m, int k)
{
    int iters = cec_c->results->iterations;
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
    PROTECT(centers_matrix = create_R_matrix(cec_c->results->centers));

    INTEGER(iterations)[0] = iters;

    for (int i = 0; i < output_size; i++)
    {
        REAL(energy_vector)[i] = cec_c->results->energy[i];
        INTEGER(clusters_number_vector)[i] = cec_c->results->clusters_number[i];
    }

    for (int i = 0; i < m; i++)
        INTEGER(assignment_vector)[i] = cec_c->results->clustering_vector[i] + 1;

    for (int i = 0; i < k; i++)
    {
        SEXP covariance;
        PROTECT(covariance = create_R_matrix(cec_c->results->covriances->mats[i]));
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
