#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"
extern "C" {
#include <R_ext/Lapack.h>
};

namespace cec {
    class utils {
    public:
        static double cec_cov_cholesky_det(mat &cov, mat &temp_mat) {
            if (cov.n == 2)
                return cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0];
            else if (!cec_cov_cholesky(cov, temp_mat))
                return std::numeric_limits<double>::quiet_NaN();
            double prod = cec_cov_diagonal_product(temp_mat);
            return handle_cholesky_nan(prod * prod);
        }

    private:
        static bool cec_cov_cholesky(const mat &cov, mat &temp_matrix) {
            int n = cov.n;
            int info;
            temp_matrix = cov;
            F77_NAME(dpotrf)("U", &n, temp_matrix.data(), &n, &info);
            return info == 0;
        }

        static double cec_cov_diagonal_product(const mat &m) {
            double res = 1.0;
            int n = m.n;
            for (int i = 0; i < n; i++)
                res *= m[i][i];
            return res;
        }

        static constexpr double ZERO_EPSILON = 1.0e-32;

        static inline double handle_zero(double d)
        {
            if (d < ZERO_EPSILON)
                return ZERO_EPSILON;
            return d;
        }

        static inline double handle_cholesky_nan(double d)
        {
            if (std::isnan(d))
                return handle_zero(0);
            return handle_zero(d);
        }

    };
};
#endif //CEC_COV_UTILS_H
