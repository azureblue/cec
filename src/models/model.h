#ifndef CEC_MODEL_H
#define CEC_MODEL_H

#include "../vec.h"

#include <cmath>

namespace cec {

    namespace math {
        const double PI = 3.14159265358979323846;
        const double E = 2.7182818284590452354;
    }

    class model {
    public:
        const int n;

        explicit model(const int n):
                n(n) {}

        virtual double cross_entropy(const mat &cov) const noexcept = 0;

        inline double energy(mat &cov, int card, int m) const noexcept {
            double p = card / (double) m;
            return p * (-std::log(p) + cross_entropy(cov));
        }
    };
}

#endif //CEC_MODEL_H
