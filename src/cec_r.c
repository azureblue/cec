#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>

#include "alloc.h"
#include "cec_r.h"
#include "cec_r_utils.h"
#include "cec.h"

static void destroy_energy_function_contexts(
	struct energy_function_context ** energy_contexts, int k);

static void destroy_energy_functions(energy_function * energy_functions);

static struct energy_function_context ** create_energy_function_contexts(
	SEXP type, SEXP params, int k, int n, int * last_error);

static energy_function * create_energy_functions(SEXP type, int k);

static SEXP create_R_result(struct cec_context * cec_c,
	struct cec_matrix * centers, int m, int k, int n);

SEXP cec_r(SEXP x, SEXP centers, SEXP iter_max, SEXP type, SEXP card_min,
	SEXP params)
{
    int m = Rf_nrows(x);
    int k = Rf_nrows(centers);
    int n = Rf_ncols(x);
    int iteration_max = Rf_asInteger(iter_max);
    int _card_min = Rf_asInteger(card_min);

    struct cec_matrix * X = create_from_R_matrix(x);
    struct cec_matrix * C = create_from_R_matrix(centers);

    if (X == NULL || C == NULL)
    {
	cec_matrix_destroy(X);
	cec_matrix_destroy(C);
	error(MALLOC_ERROR_MSG);
    }

    int last_error = 0;

    struct energy_function_context ** energy_contexts =
	    create_energy_function_contexts(type, params, k, n, &last_error);

    energy_function * energy_functions = create_energy_functions(type, k);

    struct cec_context * cec_c = create_cec_context(X, C, energy_contexts,
	    energy_functions, iteration_max, _card_min);

    if (energy_contexts == NULL || energy_functions == NULL || cec_c == NULL)
    {
	destroy_energy_function_contexts(energy_contexts, k);
	destroy_energy_functions(energy_functions);
	destroy_cec_context_results(cec_c);
	cec_matrix_destroy(X);
	cec_matrix_destroy(C);
	error(MALLOC_ERROR_MSG);
    }

    int res = cec(cec_c);

    destroy_energy_function_contexts(energy_contexts, k);
    destroy_energy_functions(energy_functions);

    SEXP result = NULL;
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
	    case POSITIVE_DEFINITE_ERROR:
		error(POSITIVE_DEFINITE_ERROR_MSG);
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

static void destroy_energy_function_contexts(
	struct energy_function_context ** energy_contexts, int k)
{
    if (energy_contexts == NULL)
	return;

    for (int i = 0; i < k; i++)
	destroy_energy_function_context(energy_contexts[i]);

    m_free(energy_contexts);

}

static void destroy_energy_functions(energy_function * energy_functions)
{
    m_free(energy_functions);
}

static energy_function * create_energy_functions(SEXP type, int k)
{
    energy_function * energy_functions = m_alloc(sizeof (energy_function) * k);
    if (energy_functions == NULL)
	return NULL;

    for (int i = 0; i < k; i++)
    {
	enum density_family family = INTEGER(type)[i];
	energy_functions[i] = energy_function_for(family);
    }

    return energy_functions;
}

static struct energy_function_context ** create_energy_function_contexts(
	SEXP type, SEXP params, int k, int n, int * last_error)
{
    struct energy_function_context ** energy_contexts = m_alloc(
	    sizeof (struct energy_function_context *) * k);

    if (energy_contexts == NULL)
	return NULL;

    for (int i = 0; i < k; i++)
    {
	SEXP param = VECTOR_ELT(params, i);
	enum density_family family = INTEGER(type)[i];
	energy_contexts[i] = create_energy_function_context(family, n);
	
	if (energy_contexts[i] == NULL)
	{
	    destroy_energy_function_contexts(energy_contexts, i + 1);
	    return NULL;
	}
	energy_contexts[i]->last_error = last_error;

	switch (family)
	{
	    case FIXED_R:
	    {
		struct context_r * c_r =
			(struct context_r *) (energy_contexts[i]->custom_context);
		c_r->r = asReal(param);
		break;
	    }
	    case GIVEN_COVARIANCE:
	    {
		struct context_gc * c_gc =
			(struct context_gc *) energy_contexts[i]->custom_context;
		copy_from_R_matrix(VECTOR_ELT(param, 0), c_gc->given_cov);
		copy_from_R_matrix(VECTOR_ELT(param, 1), c_gc->i_given_cov);		
		break;
	    }
	    case FIXEDEIGENVALUES:
	    {
		struct context_fe * c_fe =
			(struct context_fe *) energy_contexts[i]->custom_context;

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

    return energy_contexts;
}

