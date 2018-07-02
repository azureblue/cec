#ifndef ALL_H
#define ALL_H

#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class all : public model {
    public:
        explicit all(int n)
                : tmp(n, n),
                  ce_constant(n * std::log(2.0 * constants::PI * constants::E)) {}

        double cross_entropy(const mat &cov) const noexcept override {
            double det = determinant(cov, tmp);
            return (ce_constant + std::log(det)) / 2;
        }

    private:
        mutable mat tmp;
        const double ce_constant;
    };
}
#endif /* ALL_H */

