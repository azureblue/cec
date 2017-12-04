#ifndef CEC_COV_UTILS_H
#define CEC_COV_UTILS_H

#include "../vec.h"

namespace cec {
    double cec_cov_diagonal_product(const mat &m);

    double det_cholesky(mat &cov, mat &temp_mat);
}
#endif //CEC_COV_UTILS_H
