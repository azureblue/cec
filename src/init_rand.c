#include "init_rand.h"
#include "rand.h"

res_code cec_init_centers_random(const cec_mat *x, cec_mat *c) {
    int k = c->m;
    int n = c->n;
    cec_rand_init();
    for (int i = 0; i < k; i++)
    {
        double r = cec_rand();
        int p = (int) (r * x->m);
        array_copy(cec_matrix_const_row(x, p), cec_matrix_row(c, i), n);
    }
    cec_rand_end();
    return NO_ERROR;
}
