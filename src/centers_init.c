#include <string.h>
#include "centers_init.h"
#include "kmeanspp.h"
#include "init_rand.h"

res_code cec_init_centers(const cec_mat *x, cec_mat *c, init_method method) {
    switch (method) {
        case KMEANSPP:
            return cec_init_centers_kmeanspp(x, c);
        case RANDOM:
            return cec_init_centers_random(x, c);
        case NONE:
            return NO_ERROR;

        default:
            return INVALID_CENTERS_INIT_METHOD_ERROR;
    }
}

bool parse_init_method(const char * method, init_method * result) {
    if (!strcmp(method, "none"))
        return *result = NONE, true;
    if (!strcmp(method, "kmeanspp"))
        return *result = KMEANSPP, true;
    if (!strcmp(method, "random"))
        return *result = RANDOM, true;
    return false;
}
