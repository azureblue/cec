#include <string.h>
#include "centers_init.h"
#include "kmeanspp.h"
#include "init_rand.h"
#include "initial_assignment.h"

static void copy_init(initial_assignement_ctx ctx, const cec_mat *x, cec_mat *c) {
    for (int i = 0; i < c->m; i++)
        cec_matrix_copy_row(((cec_mat *) ctx), i, c, i);
}

centers_initializer *create_copy_initializer(cec_mat *centers) {
    centers_initializer * ci = alloc(centers_initializer);
    ci->ctx = cec_matrix_create_copy(centers);
    ci->init = copy_init;
    return ci;
}

void cec_init_centers(centers_initializer * ci, const cec_mat * x, cec_mat * c) {
    ci->init(ci->ctx, x, c);
}

centers_initializer * create_centers_init(cec_centers_par *centers_par, int m_max) {
    switch (centers_par->init_m) {
        case KMEANSPP:
            return create_kmeanspp_initializer(m_max);
        case RANDOM:
            return create_random_initializer();
        case NONE:
            return create_copy_initializer(centers_par->centers_mat);
    }
    return NULL;
}
