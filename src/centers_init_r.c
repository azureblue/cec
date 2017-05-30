#include "alloc.h"
#include "init_utils_r.h"
#include "cec_r_utils.h"
#include "centers_init_r.h"
#include "error_r.h"
#include "centers_init.h"
#include "rand.h"

SEXP cec_init_centers_r(SEXP x_r, SEXP k_r, SEXP method_r) {
    cec_init_env();
    cec_mat *x = create_from_R_matrix(x_r);
    cec_mat *c = cec_matrix_create(asInteger(k_r), x->n);
    init_method method;
    res_code r;
    if (parse_init_method(CHAR(STRING_ELT(method_r, 0)), &method) != true)
        r = INVALID_CENTERS_INIT_METHOD_ERROR;
    else
        r = cec_init_centers(create_from_R_matrix(x_r), c, method);

    if (r == NO_ERROR)
    {
        SEXP result;
        PROTECT(result = create_R_matrix(c));
        UNPROTECT(1);
        cec_clean_env();
        return result;
    }
    cec_clean_env();
    error_r(r);
}
