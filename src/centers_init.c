#include "centers_init.h"
#include "kmeanspp.h"
#include "init_rand.h"

cec_res cec_init_centers(const cec_mat *x, cec_mat *c, enum centers_init_method method) {
    switch (method) {
        case KMEANSPP:
            return cec_init_centers_kmeanspp(x, c);
        case RANDOM:
            return cec_init_centers_random(x, c);

        default:
            return INVALID_CENTERS_INIT_METHOD;
    }
}
