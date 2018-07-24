#ifndef DIAGONAL_H
#define DIAGONAL_H

#include "cov_utils.h"
#include "model.h"

namespace cec {
    class diagonal: public model {
    public:
        explicit diagonal(int n)
                : ce_constant(n * std::log(2.0 * m::PI * m::E)) {}

        double cross_entropy(const covariance &cov) const noexcept override {
            double diag = diagonal_product(cov);
            return (ce_constant + std::log(diag)) / 2;
        }

    private:
        const double ce_constant;
    };
}
#endif /* DIAGONAL_H */

