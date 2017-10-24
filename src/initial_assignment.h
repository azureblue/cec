#ifndef INITIAL_ASSIGNMENT_H
#define INITIAL_ASSIGNMENT_H

#include "alloc.h"
#include "matrix.h"

typedef memptr_t initial_assignement_ctx;
typedef void (*initial_assignmenet_function)(initial_assignement_ctx, const cec_mat *, const cec_mat *, vec_i *);

typedef struct {
    initial_assignmenet_function ia;
    initial_assignement_ctx ctx;
} initial_assignment;

void cec_assign_points(initial_assignment *ia, const cec_mat *x, const cec_mat *c, vec_i *result_assignment);

initial_assignment * create_initial_assignment_fixed(const vec_i *fixed_assignment);
initial_assignment * create_initial_assignment_closest();

#endif //INITIAL_ASSIGNMENT_H
