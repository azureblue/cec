#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

#include "alloc.h"
#include "cec_r.h"
#include "cec_r_utils.h"
#include "cec.h"

static void destroy_cross_entropy_contexts(
	struct cross_entropy_context ** energy_contexts, int k);

static void destroy_cross_entropy_functions(cross_entropy_function * cross_entropy_functions);

static struct cross_entropy_context ** create_cross_entropy_contexts(
	SEXP type, SEXP params, int k, int n, int * last_error);

static cross_entropy_function * create_cross_entropy_functions(SEXP type, int k);

static SEXP create_R_result(struct cec_context * cec_c,
	struct cec_matrix * centers, int m, int k, int n);

SEXP cec_r(SEXP x, SEXP centers, SEXP iter_max, SEXP type, SEXP card_min,
	SEXP params)
{
    int m = Rf_nrows(x);
    int k = Rf_nrows(centers);
    int n = Rf_ncols(x);
    int iteration_max = Rf_asInteger(iter_max);
    int card_min_int  = Rf_asInteger(card_min);

    struct cec_matrix * X = create_from_R_matrix(x);
    struct cec_matrix * C = create_from_R_matrix(centers);

    if (X == NULL || C == NULL)
    {
	cec_matrix_destroy(X);
	cec_matrix_destroy(C);
	error(MALLOC_ERROR_MSG);
    }

    int last_error = 0;

    struct cross_entropy_context ** cross_entropy_contexts =
	    create_cross_entropy_contexts(type, params, k, n, &last_error);

    cross_entropy_function * cross_entropy_functions = create_cross_entropy_functions(type, k);

    struct cec_context * cec_c = create_cec_context(X, C, cross_entropy_contexts,
	    cross_entropy_functions, iteration_max, card_min_int);

    if (cross_entropy_contexts == NULL || cross_entropy_functions == NULL || cec_c == NULL)
    {
	destroy_cross_entropy_contexts(cross_entropy_contexts, k);
	destroy_cross_entropy_functions(cross_entropy_functions);
	destroy_cec_context_results(cec_c);
	cec_matrix_destroy(X);
	cec_matrix_destroy(C);
	error(MALLOC_ERROR_MSG);
    }

    /*
     * Perform the CEC algorithm.
     */
    int res = cec(cec_c);

    destroy_cross_entropy_contexts(cross_entropy_contexts, k);
    destroy_cross_entropy_functions(cross_entropy_functions);

    SEXP result = NULL;
    
    /*
     * Prepare the results for R.
     */
    if (res == NO_ERROR)
	PROTECT(result = create_R_result(cec_c, C, m, k, n));

    destroy_cec_context_results(cec_c);

    cec_matrix_destroy(X);
    cec_matrix_destroy(C);

    if (res == NO_ERROR)
    {
	UNPROTECT(1);
	return result;
    } else
    {
	switch (res)
	{
	    case INVALID_COVARIANCE_ERROR:
		error(INVALID_COVARIANCE_ERROR_MSG);
	    case ALL_CLUSTERS_REMOVED_ERROR:
		error(ALL_CLUSTERS_REMOVED_MSG);
	    default:
		error(UNKNOWN_ERROR_MSG);
	}
    }
    return NULL;
}

static SEXP create_R_result(struct cec_context * cec_c,
	struct cec_matrix * centers, int m, int k, int n)
{
    int iters = cec_c->iterations;
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
    PROTECT(centers_matrix = create_R_matrix(centers));

    INTEGER(iterations)[0] = iters;

    for (int i = 0; i < output_size; i++)
    {
	REAL(energy_vector)[i] = cec_c->energy[i];
	INTEGER(clusters_number_vector)[i] = cec_c->clusters_number[i];
    }

    for (int i = 0; i < m; i++)
    {
	INTEGER(assignment_vector)[i] = cec_c->clustering_vector[i] + 1;
    }

    for (int i = 0; i < k; i++)
    {
	SEXP covariance;
	PROTECT(covariance = create_R_matrix(cec_c->covriances[i]));
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

static void destroy_cross_entropy_contexts(
	struct cross_entropy_context ** energy_contexts, int k)
{
    if (energy_contexts == NULL)
	return;

    for (int i = 0; i < k; i++)
	destroy_cross_entropy_context(energy_contexts[i]);

    m_free(energy_contexts);

}

static void destroy_cross_entropy_functions(cross_entropy_function * cross_entropy_functions)
{
    m_free(cross_entropy_functions);
}

static cross_entropy_function * create_cross_entropy_functions(SEXP type, int k)
{
    cross_entropy_function * cross_entropy_functions = m_alloc(sizeof (cross_entropy_function) * k);
    if (cross_entropy_functions == NULL)
	return NULL;

    for (int i = 0; i < k; i++)
    {
	enum density_family family = INTEGER(type)[i];
	cross_entropy_functions[i] = cross_entropy_for(family);
    }

    return cross_entropy_functions;
}

static struct cross_entropy_context ** create_cross_entropy_contexts(
	SEXP type, SEXP params, int k, int n, int * last_error)
{
    struct cross_entropy_context ** cross_entropy_contexts = m_alloc(
	    sizeof (struct cross_entropy_context *) * k);

    if (cross_entropy_contexts == NULL)
	return NULL;

    for (int i = 0; i < k; i++)
    {
	SEXP param = VECTOR_ELT(params, i);
	enum density_family family = INTEGER(type)[i];
	cross_entropy_contexts[i] = create_cross_entropy_context(family, n);
	
	if (cross_entropy_contexts[i] == NULL)
	{
	    destroy_cross_entropy_contexts(cross_entropy_contexts, i + 1);
	    return NULL;
	}
	cross_entropy_contexts[i]->last_error = last_error;

	switch (family)
	{
	    case FIXED_R:
	    {
		struct context_r * c_r =
			(struct context_r *) (cross_entropy_contexts[i]->custom_context);
		c_r->r = asReal(param);
		break;
	    }
	    case GIVEN_COVARIANCE:
	    {
		struct context_gc * c_gc =
			(struct context_gc *) cross_entropy_contexts[i]->custom_context;
		copy_from_R_matrix(VECTOR_ELT(param, 0), c_gc->given_cov);
		copy_from_R_matrix(VECTOR_ELT(param, 1), c_gc->i_given_cov);		
		break;
	    }
	    case FIXEDEIGENVALUES:
	    {
		struct context_fe * c_fe =
			(struct context_fe *) cross_entropy_contexts[i]->custom_context;

		array_copy(REAL(param), c_fe->given_evals, n);
		double e_prod = 1;
		for (int i = 0; i < n; i++)
		{
		    e_prod *= c_fe->given_evals[i];
		}
		c_fe->given_evals_product = e_prod;

	    }
	    case SPHERICAL:
	    case ALL:
	    case DIAGONAL:
		break;
	}
    }
    return cross_entropy_contexts;
}

