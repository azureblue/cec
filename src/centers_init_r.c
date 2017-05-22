#include "alloc.h"
#include "alloc_utils_r.h"
#include "cec_r_utils.h"
#include "centers_init_r.h"
#include "error_r.h"
#include "centers_init.h"

SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r) {
    init_mem_mg(error_r_mem_error);
    cec_mat *x = create_from_R_matrix(x_r);
    cec_mat *c = cec_matrix_create(asInteger(k_r), x->n);
    init_method method = (init_method) asInteger(method_r);
    cec_res r = cec_init_centers(create_from_R_matrix(x_r), c, method);

    if (r == NO_ERROR)
    {
        SEXP result;
        PROTECT(result = create_R_matrix(c));
        UNPROTECT(1);
        free_mem_mg();
        return result;
    }
    else {
        free_mem_mg();
        error_r(r);
    }
}
