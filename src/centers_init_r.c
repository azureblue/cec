#include "alloc.h"
#include "alloc_utils_r.h"
#include "cec_r_utils.h"
#include "centers_init_r.h"
#include "error_r.h"
#include "kmeanspp.h"
#include "rand.h"

SEXP init_kmeanspp_r(SEXP x, SEXP k)
{
    init_mem_mg(error_r_mem_error);
    struct cec_matrix * ma_x = create_from_R_matrix(x);
    struct cec_matrix * ma_c = cec_matrix_create(asInteger(k), ma_x->n);

    int r = kmeanspp(ma_x, ma_c);

    if (r == NO_ERROR)
    {
        SEXP result;
        PROTECT(result = create_R_matrix(ma_c));
        UNPROTECT(1);
        free_mem_mg();
        return result;
    } 
    else {
        free_mem_mg();
        error_r(r);
    }
}

SEXP init_random_r(SEXP x, SEXP rk)
{
    init_mem_mg(error_r_mem_error);
    int k = asInteger(rk);

    struct cec_matrix * ma_x = create_from_R_matrix(x);
    struct cec_matrix * ma_c = cec_matrix_create(k, ma_x->n);

    cec_rand_init();

    for (int i = 0; i < k; i++)
    {
        double r = cec_rand();
        int p = (int) (r * ma_x->m);
        array_copy(cec_matrix_const_row(ma_x, p), cec_matrix_row(ma_c, i), ma_x->n);
    }

    cec_rand_end();

    SEXP result;
    PROTECT(result = create_R_matrix(ma_c));
    UNPROTECT(1);
    free_mem_mg();
    return result;
}
