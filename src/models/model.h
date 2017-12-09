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
        virtual double cross_entropy(mat &cov) const = 0;

        inline double energy(mat &cov, int card, int m) const {
            double p = (card / (double) m);
            return p * (-std::log(p) + cross_entropy(cov));
        }

    };
}

#endif //CEC_MODEL_H
