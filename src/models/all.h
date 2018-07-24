#ifndef ALL_H
#define ALL_H

#include "cov_utils.h"
#include "model.h"
#include "../m.h"

namespace cec {
    class all: public model {
    public:
        explicit all(int n)
                : det_calc(n),
                  ce_constant(n * std::log(2.0 * m::PI * m::E)) {}

        double cross_entropy(const covariance &cov) const noexcept override {
            double det = det_calc.determinant(cov);
            return (ce_constant + m::log(det)) / 2;
        }

    private:
        determinant_calculator det_calc;
        const double ce_constant;
    };
}
#endif /* ALL_H */

