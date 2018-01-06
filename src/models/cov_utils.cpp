#include "cov_utils.h"
#include "constants.h"

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

static bool cholesky(const cec::mat &cov, cec::mat &tmp) {
    int n = cov.n;
    int info;
    tmp = cov;
    dpotrf_("U", &n, tmp.data(), &n, &info);
    return info == 0;
}

double cec::diagonal_product(const mat &cov) {
    int n = cov.n;
    double res = 1.0;
    for (int i = 0; i < n; i++)
        res *= cov[i][i];
    return res;
}

double cec::trace(const mat &cov) {
    int n = cov.n;
    double tr = 0;
    for (int i = 0; i < n; i++)
        tr += cov[i][i];
    return tr;
}

double cec::determinant(const mat &cov, mat &tmp_mat) {
    if (cov.n == 1)
        return cov[0][0];
    if (cov.n == 2)
        return cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
    if (!cholesky(cov, tmp_mat))
        return constants::QNAN;
    double prod = diagonal_product(tmp_mat);
    return handle_cholesky_nan(prod * prod);
}
