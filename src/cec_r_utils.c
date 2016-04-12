#include "cec_r_utils.h"

struct cec_matrix * create_from_R_matrix(SEXP R_ma)
{
    int m = Rf_nrows(R_ma);
    int n = Rf_ncols(R_ma);

    struct cec_matrix * ma = cec_matrix_create(m, n);

    if (!ma)
        return NULL;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            cec_matrix_set_element(ma, i, j, REAL(R_ma)[j * m + i]);

    return ma;
}

void copy_from_R_matrix(SEXP R_ma, struct cec_matrix * ma)
{
    int m = Rf_nrows(R_ma);
    int n = Rf_ncols(R_ma);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            cec_matrix_set_element(ma, i, j, REAL(R_ma)[j * m + i]);
}

SEXP create_R_matrix(struct cec_matrix * ma)
{
    int m = ma->m;
    int n = ma->n;

    SEXP R_ma;

    PROTECT(R_ma = allocMatrix(REALSXP, m, n));

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            REAL(R_ma)[j * m + i] = cec_matrix_element(ma, i, j);

    UNPROTECT(1);

    return R_ma;
}
