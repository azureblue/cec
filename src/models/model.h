#ifndef CEC_MODEL_H
#define CEC_MODEL_H

#include "../vec.h"

namespace cec {

    class model {
    public:
        virtual double cross_entropy(mat &cov) const = 0;

        inline double energy(mat &cov, int card, int m) const {
            double p = (card / (double) m);
            return p * (-log(p) + cross_entropy(cov));
        }

    };
}

#endif //CEC_MODEL_H
