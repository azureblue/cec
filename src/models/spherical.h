#ifndef SPHERICAL_H
#define SPHERICAL_H

#include "cov_utils.h"
#include "model.h"

namespace cec {
    class spherical: public model {
    public:
        explicit spherical(int n)
                : n(n),
                  ce_constant(std::log(2.0 * m::PI * m::E / n)) {}

        double cross_entropy(const covariance &cov) const noexcept override {
            double tr = trace(cov);
            return (ce_constant + std::log(tr)) * n / 2;
        }

    private:
        const int n;
        const double ce_constant;
    };
}
#endif /* SPHERICAL_H */

