#ifndef MEAN_H
#define MEAN_H

#include "cov_utils.h"
#include "model.h"
#include "../m.h"

namespace cec {
    class fixed_mean: public model {
    public:
        explicit fixed_mean(int n, const row &fixed_mean)
                : det_calc(n),
                  mahalanobis_dist_calc(n),
                  cov_inv(n, n),
                  mean(fixed_mean),
                  ce_constant(n * std::log(2.0 * m::PI * m::E)) {}

        double cross_entropy(const covariance &cov) const noexcept override {
            double det = det_calc.determinant(cov);
            if (!invert(cov, cov_inv))
                return m::QNAN;
            double md2 = mahalanobis_dist_calc.mahalanobis2(cov_inv, cov.mean(), mean);
            return (m::log(1 + md2) + ce_constant + m::log(det)) / 2;
        }

    private:
        determinant_calculator det_calc;
        mahalanobis_dist_calculator mahalanobis_dist_calc;
        mutable mat cov_inv;
        const vec mean;
        const double ce_constant;
    };
}
#endif /* MEAN_H */
