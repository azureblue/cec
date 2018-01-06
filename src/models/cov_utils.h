#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"

namespace cec {
    double diagonal_product(const mat &cov);

    double determinant(const mat &cov, mat &tmp_mat);

    double trace(const mat &cov);
}
#endif //CEC_COV_UTILS_H
