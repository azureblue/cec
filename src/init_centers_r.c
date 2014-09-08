#include "init_centers_r.h"
#include "cec_r_utils.h"
#include "kmeanspp.h"
#include "errors.h"
#include <stdio.h>

SEXP init_kmeanspp_r(SEXP x, SEXP k)
{
    struct cec_matrix * ma_x = create_from_R_matrix(x);

    if (ma_x == NULL)
	error(MALLOC_ERROR_MSG);

    struct cec_matrix * ma_c = cec_matrix_create(asInteger(k), ma_x->n);

    if (ma_c == NULL)
    {
	cec_matrix_destroy(ma_x);
	error(MALLOC_ERROR_MSG);
    }

    int r = kmeanspp(ma_x, ma_c);

    if (r == NO_ERROR)
    {
	SEXP result;
	PROTECT(result = create_R_matrix(ma_c));
	cec_matrix_destroy(ma_x);
	cec_matrix_destroy(ma_c);
	UNPROTECT(1);
	return result;
    } else
    {
	cec_matrix_destroy(ma_x);
	cec_matrix_destroy(ma_c);
	switch (r)
	{
	    case MALLOC_ERROR:
		error(MALLOC_ERROR_MSG);
	    default:
		error(UNKNOWN_ERROR_MSG);
	}
    }
}

SEXP init_random_r(SEXP x, SEXP rk)
{
    int k = asInteger(rk);

    struct cec_matrix * ma_x = create_from_R_matrix(x);

    if (ma_x == NULL)
	error(MALLOC_ERROR_MSG);

    struct cec_matrix * ma_c = cec_matrix_create(k, ma_x->n);

    if (ma_c == NULL)
    {
	cec_matrix_destroy(ma_x);
	error(MALLOC_ERROR_MSG);
    }

    for (int i = 0; i < k; i++)
    {
	double r = (double) rand();
	int p = (int) ((r / (double) RAND_MAX) * ma_x->m);
	array_copy(cec_matrix_row(ma_x, p), cec_matrix_row(ma_c, i), ma_x->n);
    }

    SEXP result;
    PROTECT(result = create_R_matrix(ma_c));
    cec_matrix_destroy(ma_x);
    cec_matrix_destroy(ma_c);
    UNPROTECT(1);
    return result;
}
