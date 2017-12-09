#include "cov_utils.h"

extern "C" {
#include <R_ext/Lapack.h>
}

const double ZERO_EPSILON = 1.0e-32;

static inline double handle_zero(double d) {
    if (d < ZERO_EPSILON)
        return ZERO_EPSILON;
    return d;
}

static inline double handle_cholesky_nan(double d) {
    if (std::isnan(d))
        return handle_zero(0);
    return handle_zero(d);
}

static bool cec_cov_cholesky(const cec::mat &cov, cec::mat &temp_matrix) {
    int n = cov.n;
    int info;
    temp_matrix = cov;
    F77_NAME(dpotrf)("U", &n, temp_matrix.data(), &n, &info);
    return info == 0;
}

double cec::cec_cov_diagonal_product(const cec::mat &m) {
    double res = 1.0;
    int n = m.n;
    for (int i = 0; i < n; i++)
        res *= m[i][i];
    return res;
}

double cec::det_cholesky(cec::mat &cov, cec::mat &temp_mat) {
    if (cov.n == 2)
        return cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
    else if (!cec_cov_cholesky(cov, temp_mat))
        return std::numeric_limits<double>::quiet_NaN();
    double prod = cec_cov_diagonal_product(temp_mat);
    return handle_cholesky_nan(prod * prod);
}
