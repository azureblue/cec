#include "cov_utils.h"
#include "../m.h"

#include <R_ext/Lapack.h>

const double ZERO_EPSILON = 1.0e-32;

static inline double handle_zero(double d) {
    if (d < ZERO_EPSILON)
        return ZERO_EPSILON;
    return d;
}

static inline double handle_cholesky_nan(double d) {
    if (cec::m::isnan(d))
        return handle_zero(0);
    return handle_zero(d);
}

static bool cholesky(const cec::mat &cov, cec::mat &dst) {
    int n = cov.n;
    int info;
    dst = cov;
    dpotrf_("U", &n, dst.data(), &n, &info);
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

void cec::multiply(const mat &a, const mat &b, mat &dst) {
    int n = a.n;
    if (n < 8) {
        dst.fill(0.0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++)
                    dst[i][j] += a[k][i] * b[j][k];
        return;
    }
    double zero = 0;
    double one = 1;
    dsymm_("L", "L", &n, &n, &one, a.data(), &n, b.data(), &n, &zero, dst.data(), &n);
}

bool cec::invert(const mat &cov, mat &dst) {
    int n = cov.n;
    int info;
    if (!cholesky(cov, dst))
        return false;
    dpotri_("U", &n, dst.data(), &n, &info);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            dst[i][j] = dst[j][i];
    return info == 0;
}

bool cec::eigenvalues_calculator::eigenvalues(const cec::mat &cov, double *res) const noexcept {
    int n = cov.n;
    int info;
    tmp = cov;
    dsyev_("N", "U", &n, tmp.data(), &n, res, workspace.data(), &workspace.size, &info);
    return info == 0;
}

double cec::determinant_calculator::determinant(const cec::mat &cov) const noexcept {
    if (cov.n == 1)
        return cov[0][0];
    if (cov.n == 2)
        return cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
    if (!cholesky(cov, tmp))
        return m::QNAN;
    double prod = diagonal_product(tmp);
    return handle_cholesky_nan(prod * prod);
}

double
cec::mahalanobis_dist_calculator::mahalanobis2(const mat &cov_inv, const row &mean, const row &x) const {
    int n = cov_inv.n;
    for (int i = 0; i < n; i++)
        tmp[i] = x[i] - mean[i];
    double res = 0;
    for (int i = 0; i < n; i++) {
        double acc = 0.0;
        for (int j = 0; j < n; j++)
            acc += tmp[j] * cov_inv[j][i];
        res += acc * tmp[i];
    }
    return res;
}
