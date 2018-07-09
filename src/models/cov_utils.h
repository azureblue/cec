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
                : n(n),
                  tmp(n, n),
                  workspace(WORKSPACE_SIZE_MUL * n) {}

        bool eigenvalues(const mat &cov, double *res) const noexcept;

    private:
        const int n;
        mutable mat tmp;
        mutable std::vector<double> workspace;
        static const int WORKSPACE_SIZE_MUL = 130;
    };

}
#endif //CEC_COV_UTILS_H
