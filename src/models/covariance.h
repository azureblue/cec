#ifndef COVARIANCE_H
#define COVARIANCE_H

#include "cov_utils.h"
#include "model.h"
#include "../exceptions.h"

namespace cec {
    class covariance : public model {
    public:
        explicit covariance(int n, mat cov)
                : cov_inv(std::move(inv(cov))),
                  tmp(n, n),
                  ce_constant(std::log(std::pow(2.0 * m::PI, n) * det(cov)) / 2.0) {}

        double cross_entropy(const mat &cov) const noexcept override {
            multiply(cov_inv, cov, tmp);
            double tr = trace(tmp);
            return ce_constant + tr / 2;
        }

    private:
        const mat cov_inv;
        mutable mat tmp;
        const double ce_constant;
        
        static mat inv(const mat& cov) {
            mat dst(cov);
            if (!invert(cov, dst))
                throw new invalid_model_parameter("invalid covariance (not positive definite)");
            return dst;
        }

        static double det(const mat& cov) {
            mat tmp(cov);
            double det = determinant(cov, tmp);
            if (std::isnan(det))
                throw new invalid_model_parameter("invalid covariance (not positive definite)");
            return det;
        }
    };
}
#endif /* COVARIANCE_H */

