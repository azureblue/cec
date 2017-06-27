#ifndef CEC_TYPE_H
#define CEC_TYPE_H

#include <stdbool.h>

enum density_family {
    ALL, COVARIANCE, DIAGONAL, EIGENVALUES, FIXED_R, SPHERICAL
};

bool cec_parse_type(const char * type, enum density_family * result);

#endif //CEC_TYPE_H
