#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include <vector>

#include "../vec.h"

namespace cec {

    double diagonal_product(const mat &cov);

    double determinant(const mat &cov, mat &tmp);

    bool invert(const mat &cov, mat &dst);

    void multiply(const mat &a, const mat &b, mat &dst);

    double trace(const mat &cov);

    class eigenvalues_calculator {
    public:
        explicit eigenvalues_calculator(const int n)
                : tmp(n, n),
                  workspace(WORKSPACE_SIZE_MUL * n) {}

        bool eigenvalues(const mat &cov, double *res) const noexcept;

    private:
        mutable mat tmp;
        mutable vec workspace;
        static const int WORKSPACE_SIZE_MUL = 130;
    };

    class determinant_calculator {
    public:
        explicit determinant_calculator(const int n)
                : tmp(n, n) {}

        double determinant(const mat &cov) const noexcept;

    private:
        mutable mat tmp;
    };

    class mahalanobis_dist_calculator {
    public:
        explicit mahalanobis_dist_calculator(const int n)
                : tmp(n) {}

        double mahalanobis2(const mat &cov_inv, const row &mean, const row &x) const;

    private:
        mutable vec tmp;
    };
}
#endif //CEC_COV_UTILS_H
