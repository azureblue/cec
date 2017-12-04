#ifndef ALL_H
#define ALL_H
#include "../vec.h"
#include "cov_utils.h"
#include "model.h"

namespace cec {

    class all: public model {
    public:
        explicit all(int n) : temp_mat(n, n) {}

        double cross_entropy(mat &cov) const {
            double det = det_cholesky(cov, const_cast<mat &>(temp_mat));
            return (cov.n / 2.0) * log(2.0 * M_PI * M_E) + (1.0 / 2.0) * log(det);
        }

    private:
        mat temp_mat;
    };
}
#endif /* ALL_H */

