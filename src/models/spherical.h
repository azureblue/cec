#ifndef SPHERICAL_H
#define SPHERICAL_H

#include "../vec.h"
#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class spherical : public model {
    public:
        explicit spherical(int n)
                : n(n),
                  spherical_const(std::log(2.0 * constants::PI * constants::E / n)) {}

        double cross_entropy(const mat &cov) const noexcept {
            double tr = trace(cov);
            return (spherical_const + std::log(tr)) * n / 2;
        }

    private:
        const int n;
        const double spherical_const;
    };
}
#endif /* SPHERICAL_H */

