#ifndef FIXED_RADIUS_H
#define FIXED_RADIUS_H

#include "cov_utils.h"
#include "model.h"
#include "constants.h"

namespace cec {
    class fixed_radius : public model {
    public:
        explicit fixed_radius(int n, double r)
                : r(r),
                  ce_const(n * std::log(2.0 * constants::PI * r) / 2.0) {}

        double cross_entropy(const mat &cov) const noexcept {
            double tr = trace(cov);
            return ce_const + tr / (2.0 * r);
        }

    private:
        const double r;
        const double ce_const;
    };
}
#endif /* FIXED_RADIUS_H */

