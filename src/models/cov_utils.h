#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"

namespace cec {
    double diagonal_product(const mat &cov);

    double determinant(const mat &cov, mat &tmp);

    bool invert(const mat &cov, mat &dst);

    void multiply(const mat &m1, const mat &m2, mat &dst);

    double trace(const mat &cov);

}
#endif //CEC_COV_UTILS_H
