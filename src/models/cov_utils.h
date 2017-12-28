#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"

namespace cec {
    double diagonal_product(const mat &m);

    double det_cholesky(const mat &cov, mat &tmp_mat);
}
#endif //CEC_COV_UTILS_H
