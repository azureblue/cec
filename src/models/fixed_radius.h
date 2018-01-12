#ifndef FIXED_RADIUS_H
#define FIXED_RADIUS_H

#include "cov_utils.h"
#include "model.h"

namespace cec {
    class fixed_radius : public model {
    public:
        explicit fixed_radius(int n, double r)
                : r(r),
                  ce_constant(n * m::log(2.0 * m::PI * r) / 2.0) {}

        double cross_entropy(const mat &cov) const noexcept override {
            double tr = trace(cov);
            return ce_constant + tr / (2.0 * r);
        }

    private:
        const double r;
        const double ce_constant;
    };
}
#endif /* FIXED_RADIUS_H */

