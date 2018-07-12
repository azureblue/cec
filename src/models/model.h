#ifndef CEC_MODEL_H
#define CEC_MODEL_H

#include "../vec.h"
#include "../m.h"
#include "../cov.h"

namespace cec {
    class model {
    public:
        explicit model() = default;

        virtual ~model() = default;

        virtual double cross_entropy(const covariance &cov) const noexcept = 0;

        inline double energy(const covariance &cov, int m) const noexcept {
            double p = cov.card() / (double) m;
            return p * (-m::log(p) + cross_entropy(cov));
        }
    };
}

#endif //CEC_MODEL_H
