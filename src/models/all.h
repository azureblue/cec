#ifndef ALL_H
#define ALL_H

#include "../vec.h"
#include "cov_utils.h"
#include "model.h"

namespace cec {

    class all: public model {
    public:
        explicit all(int n) : model(n), temp_mat(n, n) {}

        double cross_entropy(const mat &cov) const noexcept {
            double det = det_cholesky(cov, temp_mat);
            return (n / 2.0) * std::log(2.0 * math::PI * math::E) + (1.0 / 2.0) * std::log(det);
        }

    private:
        mutable mat temp_mat;
    };
}
#endif /* ALL_H */

