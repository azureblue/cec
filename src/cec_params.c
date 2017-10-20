#include <string.h>
#include "cec_params.h"

bool parse_init_method(const char *method, enum centers_init_method *result) {
    if (!strcmp(method, "none"))
        return *result = NONE, true;
    if (!strcmp(method, "kmeanspp"))
        return *result = KMEANSPP, true;
    if (!strcmp(method, "random"))
        return *result = RANDOM, true;
    return false;
}

bool cec_parse_type(const char * type, enum density_family *result) {
    if (!strcmp(type, "all"))
        return *result = ALL, true;
    if (!strcmp(type, "covariance"))
        return *result = COVARIANCE, true;
    if (!strcmp(type, "diagonal"))
        return *result = DIAGONAL, true;
    if (!strcmp(type, "eigenvalues"))
        return *result = EIGENVALUES, true;
    if (!strcmp(type, "fixed_r") || !strcmp(type, "fixedr"))
        return *result = FIXED_R, true;
    if (!strcmp(type, "spherical"))
        return *result = SPHERICAL, true;
    return false;
}
