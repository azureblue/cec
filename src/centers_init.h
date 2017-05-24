#include <stdbool.h>
#include "matrix.h"

#ifndef CEC_CENTERS_INIT_H
#define CEC_CENTERS_INIT_H

enum centers_init_method {
    NONE = 0,
    KMEANSPP = 1,
    RANDOM = 2
};

#define UNDEFINED -1

typedef enum centers_init_method init_method;

res_code cec_init_centers(const cec_mat * x, cec_mat * c, init_method method);

bool parse_init_method(const char * method, init_method * result);

#endif //CEC_CENTERS_INIT_H
