#include "matrix.h"

#ifndef CEC_CENTERS_INIT_H
#define CEC_CENTERS_INIT_H

enum centers_init_method {
    KMEANSPP = 0,
    RANDOM = 1
};

typedef enum centers_init_method init_method;

cec_res cec_init_centers(const cec_mat * x, cec_mat * c, init_method method);

#endif //CEC_CENTERS_INIT_H
