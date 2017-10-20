#include "init_rand.h"
#include "rand.h"

static void init(centers_init_ctx ctx, const cec_mat *x, cec_mat *c) {
    int k = c->m;
    int n = c->n;
    for (int i = 0; i < k; i++) {
        double r = cec_rand();
        int p = (int) (r * x->m);
        array_copy(cec_matrix_const_row(x, p), cec_matrix_row(c, i), n);
    }
}

centers_init *create_random_initializer() {
    centers_init *ci = alloc(centers_init);
    ci->ctx = NULL;
    ci->init = init;
    return ci;
}
