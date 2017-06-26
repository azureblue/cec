#include <math.h>
#include <stdbool.h>
#include <R_ext/Lapack.h>
#include "cov_utils.h"

double cec_cov_trace(const cec_mat *m)
{
    double res = 0;
    for (int i = 0; i < m->n; i++)
        res += cec_matrix_element(m, i, i);
    return res;
}

double cec_cov_diagonal_product(const cec_mat *m)
{
    double res = 1.0;
    int n = m->n;
    for (int i = 0; i < n; i++)
        res *= cec_matrix_element(m, i, i);
    return res;
}

void cec_cov_multiply(const cec_mat *m1, const cec_mat *m2, cec_mat *restrict dest)
{
    int n = m1->n;
    cec_matrix_set(dest, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                cec_matrix_set_element(dest, i, j,
                    cec_matrix_element(dest, i, j)
                    + cec_matrix_element(m1, k, i)
                    * cec_matrix_element(m2, j, k));
}

bool cec_cov_eigenvalues(const cec_mat *sym_matrix, cec_mat *temp_matrix,
        cec_mat *workspace, double *values) {
    int n = sym_matrix->n;
    int info;
    array_copy(sym_matrix->data, temp_matrix->data, n * n);
    int workspace_len = workspace->m * workspace->n;
    F77_NAME(dsyev)("N", "U", &n, temp_matrix->data, &n, values, workspace->data, &workspace_len, &info);
    return info == 0;
}

bool cec_cov_cholesky(const cec_mat *restrict sym_matrix, cec_mat *restrict temp_matrix) {
    int n = sym_matrix->n;
    int info;
    array_copy(sym_matrix->data, temp_matrix->data, n * n);
    F77_NAME(dpotrf)("U", &n, temp_matrix->data, &n, &info);
    return info == 0;
}

double cec_cov_cholesky_det(const cec_mat *m, cec_mat *temp) {
    if (m->n == 2)
        return m->data[0] * m->data[3] - m->data[1] * m->data[2];
    else if (!cec_cov_cholesky(m, temp))
            return NAN;
    double prod = cec_cov_diagonal_product(temp);
    return prod * prod;
}

static inline void cec_cov_change(cec_mat *restrict dest_covariance, const cec_mat *restrict covariance,
        const double *restrict mean, double const *restrict point, double cov_mul, double new_cov_point_mul)
{
    const int n = covariance->n;

    for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++)
        {
            double val = cec_matrix_element(covariance, j, k) * cov_mul
                    + (mean[j] - point[j]) * (mean[k] - point[k]) * new_cov_point_mul;

            cec_matrix_set_element(dest_covariance, j, k, val);
        }
}

void cec_cov_add_point(cec_mat *restrict dest_covariance, const cec_mat *restrict covariance,
        const double *restrict mean, double const *restrict point, int card)
{
    double card_n = card + 1;
    cec_cov_change(dest_covariance, covariance, mean, point, card / card_n, card / (card_n * card_n));
}

void cec_cov_remove_point(cec_mat *restrict dest_covariance, const cec_mat * restrict covariance,
        const double *restrict mean, double const *restrict point, int card)
{
    double card_n = card - 1;
    cec_cov_change(dest_covariance, covariance, mean, point, card / card_n, -card / (card_n * card_n));
}
