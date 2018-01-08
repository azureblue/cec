#ifndef DIAGONAL_H
#define DIAGONAL_H

#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class diagonal : public model {
    public:
        explicit diagonal(int n)
                : ce_constant(n * std::log(2.0 * constants::PI * constants::E)) {}

        double cross_entropy(const mat &cov) const noexcept override {
            double diag = diagonal_product(cov);
            return (ce_constant + std::log(diag)) / 2;
        }

    private:
        const double ce_constant;
    };
}
#endif /* DIAGONAL_H */

