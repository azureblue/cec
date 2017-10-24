#include "initial_assignment.h"
#include "cec.h"

void assign_fixed(initial_assignement_ctx ctx, const cec_mat *x, const cec_mat *c, vec_i *result_assignment) {
    vec_i_copy(((vec_i *) ctx), result_assignment);
}

void assign_closest(initial_assignement_ctx ctx, const cec_mat *x, const cec_mat *c, vec_i *result_assignment) {
    int m = x->m;
    int k = c->m;
    int n = x->n;
    for (int i = 0; i < m; i++) {
        double dist = BIG_DOUBLE;
        for (int j = 0; j < k; j++) {
            double dist_temp = dist_sq(mat_crow(x, i), mat_crow(c, j), n);
            if (dist > dist_temp) {
                dist = dist_temp;
                result_assignment->ar[i] = j;
            }
        }
    }
}

initial_assignment *create_initial_assignment_closest() {
    initial_assignment *ia = alloc(initial_assignment);
    ia->ia = assign_closest;
    return ia;
}

initial_assignment *create_initial_assignment_fixed(const vec_i *fixed_assignment) {
    initial_assignment *ia = alloc(initial_assignment);
    ia->ia = assign_fixed;
    ia->ctx = vec_i_create(fixed_assignment->len);
    vec_i_copy(fixed_assignment, ((vec_i *) ia->ctx));
    return ia;
}

void cec_assign_points(initial_assignment *ia, const cec_mat *x, const cec_mat *c, vec_i *result_assignment) {
    ia->ia(ia->ctx, x, c, result_assignment);
}