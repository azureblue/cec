#ifndef DIAGONAL_H
#define DIAGONAL_H

#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class diagonal : public model {
    public:
        explicit diagonal(int n)
                : diagonal_const(n * std::log(2.0 * constants::PI * constants::E)) {}

        double cross_entropy(const mat &cov) const noexcept {
            double diag = diagonal_product(cov);
            return (diagonal_const + std::log(diag)) / 2;
        }

    private:
        const double diagonal_const;
    };
}
#endif /* DIAGONAL_H */

