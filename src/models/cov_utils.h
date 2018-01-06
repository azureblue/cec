#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"

namespace cec {

    double diagonal_product(const mat &cov);

    double determinant(const mat &cov, mat &tmp);

    bool invert(const mat &cov, mat &dst);

    void multiply(const mat &m1, const mat &m2, mat &dst);

    double trace(const mat &cov);

    class eigenvalues_calculator {
    public:
        explicit eigenvalues_calculator(const int n)
                : n(n),
                  tmp(n, n),
                  workspace(130 * n) {}

        bool eigenvalues(const mat &cov, double *res) const noexcept;

    private:
        const int n;
        mutable mat tmp;
        mutable std::vector<double> workspace;
    };

}
#endif //CEC_COV_UTILS_H
