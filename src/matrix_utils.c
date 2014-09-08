#include "matrix_utils.h"
#include <R_ext/Lapack.h>

int cec_eigenvalues_sm(const struct cec_matrix * sym_matrix, struct cec_matrix * temp_matrix, double * work, const int lwork, double * values)
{
    int n = sym_matrix->n;
    int info;
    array_copy(sym_matrix->data, temp_matrix->data, n * n);
    F77_NAME(dsyev)("N", "U", &n, temp_matrix->data, &n, values, work, &lwork, &info);
    return info;
}