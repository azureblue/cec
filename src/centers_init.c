#include <string.h>
#include "centers_init.h"
#include "kmeanspp.h"
#include "init_rand.h"

static void no_init(centers_init_ctx ctx, const cec_mat *x, cec_mat *c) {
}

static centers_init *create_none_initializer() {
    centers_init * ci = alloc(centers_init);
    ci->init = no_init;
    return ci;
}

void cec_init_centers(centers_init * ci, const cec_mat * x, cec_mat * c) {
    ci->init(ci->ctx, x, c);
}

centers_init * create_centers_init(cec_centers_par *centers_par, int m_max) {
    switch (centers_par->init_m) {
        case KMEANSPP:
            return create_kmeanspp_initializer(m_max);
        case RANDOM:
            return create_random_initializer();
        case NONE:
            return create_none_initializer();
    }
    return NULL;
}
