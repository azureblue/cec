#ifndef ALL_H
#define ALL_H

#include "../vec.h"
#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class all : public model {
    public:
        explicit all(int n)
                : tmp(n, n),
                  all_const(n * std::log(2.0 * constants::PI * constants::E)) {}

        double cross_entropy(const mat &cov) const noexcept {
            double det = determinant(cov, tmp);
            return (all_const + std::log(det)) / 2;
        }

    private:
        mutable mat tmp;
        const double all_const;
    };
}
#endif /* ALL_H */

